suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("ggbeeswarm")
  library("here")
  library("BASiCS")
  library("coda")
  library("viridis")
})

source(here("src/analysis/functions.R"))
theme_set(theme_bw())

datasets <- c("buettner", "chen", "tung", "zeisel")
data_dims <- vapply(
    datasets,
    function(x) {
        suppressMessages(
            dim(
                readRDS(paste0("rdata/", x, ".rds"))
            )
        )
    },
    FUN.VALUE = numeric(2)
)
data_dims <- as.data.frame(t(data_dims))
colnames(data_dims) <- c("nGenes", "nCells")
data_dims[["data"]] <- datasets
data_dims[["libsize"]] <- vapply(datasets,
    function(x) {
        median(colSums(counts(readRDS(paste0("rdata/", x, ".rds")))))
    },
    numeric(1)
)
