library("SingleCellExperiment")
library("argparse")
library("xtable")

parser <- ArgumentParser()
parser$add_argument("-i", "--input")

args <- parser$parse_args()

sces <- args[["input"]]
if (is.null(sces)) {
    sces <- list.files("rdata", full.names = TRUE, pattern="*.rds")
}

dims <- lapply(sces, function(x) dim(readRDS(x)))

df <- data.frame(Dataset = gsub("rdata/", "", sces))
# df$Dataset <- gsub(".rds", "", df$Dataset)
# df$Dataset <- gsub("([\\w])([\\w]+)", "\\U\\1\\L\\2", df$Dataset, perl = TRUE)
df$Dataset <- NULL
df$"Original publication" <- c("Buettner2015", "Chen2017a", "Ibarra-Soria2018", "Tung2017", "Zeisel2015")
df$"Original publication" <- paste0("\\cite{", df$"Original publication", "}")
df <- cbind(df, as.data.frame(do.call(rbind, dims)))
colnames(df)[2:3] <- c("Genes", "Cells")
df$"Cell type" <- c("T (G2M only)", "Astrocytes", "Mesoderm", "iPSC", "CA1+ pyramidal")
df$Type <- c("Plate", "Droplet", "Droplet", "Plate", "Plate")
df$UMI <- c("No", "Yes", "Yes", "Yes", "Yes")


out <- xtable(
    df,
    align = "rr|rrrrr",
    caption = "A summary of the datasets used in the testing of scalable inference methods.",
    label = "tab:CellsGenes"
)
print(out,
    include.rownames = FALSE,
    hline.after = 0,
    # latex.environments = "widestuff",
    sanitize.text.function = function(x) {x},
    file = "tables/data-summary.tex"
)
