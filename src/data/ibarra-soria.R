library("scran")
library("scater")
library("SingleCellExperiment")
options(timeout = 10000)

if (!file.exists("downloads/rawCounts.tsv")) {
    file <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6153/rawCounts.tsv"
    download.file(
        file,
        destfile = "downloads/rawCounts.tsv"
    )
}
rawCounts <- read.delim("downloads/rawCounts.tsv", header = TRUE)

if (!file.exists("downloads/cellAnnotation.tsv")) {
    file <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6153/cellAnnotation.tsv"
    download.file(
        file,
        destfile = "downloads/cellAnnotation.tsv"
    )
}

cluster_labels <- read.table("downloads/cellAnnotation.tsv",
    sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cluster_labels[["Cell_type"]] <- cluster_labels$cellType
cluster_labels[["Cell_type"]] <- sub(
    "presomiticMesoderm",
    "PSM",
    cluster_labels[["Cell_type"]]
)
cluster_labels[["Cell_type"]] <- sub("somiticMesoderm",
    "SM", 
    cluster_labels[["Cell_type"]]
)

ind_som <- which(cluster_labels[["Cell_type"]] == "PSM" |
    cluster_labels[["Cell_type"]] == "SM")
rawCounts <- rawCounts[, ind_som]
cluster_labels <- cluster_labels[ind_som, ]

droplet_sce <- SingleCellExperiment(
    assays = list(counts = as(as.matrix(rawCounts), "dgCMatrix"))
)
rm(rawCounts)
colData(droplet_sce) <- DataFrame(
    cluster_labels,
    subCellType = sub("_.*", "", cluster_labels$cell)
)
ind_expressed <- which(Matrix::rowMeans(counts(droplet_sce)) > 0.1)
droplet_sce <- droplet_sce[ind_expressed, ]

droplet_sce <- computeSumFactors(droplet_sce,
    clusters = colData(droplet_sce)$subCellType
)
droplet_sce <- logNormCounts(droplet_sce)
droplet_sce <- runPCA(droplet_sce)


# Cell types identified by clustering
# plotReducedDim(droplet_sce, dimred = "PCA", colour_by = "subCellType") +
#   scale_fill_manual(name = "cellType", values = c("coral4", "steelblue", "limegreen"))

## exclude outliers (TODO: justify)
ind_retain <- colData(droplet_sce)$subCellType != "presomiticMesoderm.b"
droplet_sce <- droplet_sce[, ind_retain]

counts(droplet_sce) <- as.matrix(counts(droplet_sce))

droplet_sce$BatchInfo <- as.character(round(droplet_sce$sample))
colnames(droplet_sce) <- paste("Cell", seq_len(ncol(droplet_sce)))

saveRDS(droplet_sce, "rdata/ibarra-soria.rds")
