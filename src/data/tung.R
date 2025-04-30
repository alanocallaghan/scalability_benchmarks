library("scater")
library("BASiCS")
library("here")

data.dir <- "downloads"
tungfile <- file.path(data.dir, "tung_molecules-filter.txt")
erccfile <- file.path(data.dir, "tung_ercc_conc.txt")
system(
    paste(
        "wget https://raw.githubusercontent.com/jdblischak/singleCellSeq/b6ebd4a36925995e0b9b806a4172e9fb1f99e7b6/data/expected-ercc-molecules.txt",
        "-O", erccfile
    )
)


system(
    paste(
        "wget", 
        "https://github.com/Alanocallaghan/singleCellSeq/raw/master/data/molecules-filter.txt",
        "-O", tungfile
    )
)

reads <- read.delim(
    tungfile,
    stringsAsFactors = FALSE
)
rownames(reads) <- reads[[1]]
reads <- reads[, -1]
reads <- as.matrix(reads)

cell_lines <- gsub("(.*)\\.r\\d\\..*", "\\1", colnames(reads))
cell_ind <- 1
ind_cell_one <- cell_lines == sort(unique(cell_lines))[[cell_ind]]

reads <- reads[, ind_cell_one]

libsize_drop <- isOutlier(colSums(reads), nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(colSums(reads != 0), nmads = 3, type = "lower", log = TRUE)

reads <- reads[, !(libsize_drop | feature_drop)]

spikes <- rownames(reads)[grep("ERCC", rownames(reads))]

ind_expressed <- rowMeans(reads) >= 1 & rowMeans(reads != 0) > 0.5
reads <- reads[ind_expressed, ]


ERCC_num <- read.table(
    erccfile,
    header = TRUE,
    sep = "\t",
    fill = TRUE
)
ERCC_num <- ERCC_num[, c("id", "ercc_molecules_well")]
rownames(ERCC_num) <- ERCC_num[["id"]]

SpikeInput <- ERCC_num[grep("ERCC", rownames(reads), value = TRUE), ]

batches <- gsub(".*\\.(r\\d)\\..*", "\\1", colnames(reads))

bd <- newBASiCS_Data(
    Counts = reads,
    BatchInfo = batches,
    Tech = rownames(reads) %in% spikes,
    SpikeInfo = SpikeInput
)
colnames(bd) <- colnames(reads)
saveRDS(bd, file = "rdata/tung.rds")
