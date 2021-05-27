library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("coda")

theme_set(theme_bw())
source(here("src/analysis/functions.R"))


source(here("src/analysis/data_comparison.R"))


# advi_files <- list.files("outputs/advi", full.names = TRUE)
# advi_triplets <- file2triplets(advi_files)
# advi_elbo <- lapply(advi_triplets, function(x) readRDS(x[[3]]))
# advi_triplets <- lapply(advi_triplets, function(x) x[-3])
# advi_df <- read_triplets(advi_triplets)


# source(here("src/analysis/elbo_plots.R"))

dc_files <- list.files("outputs/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(file2triplets(dc_files), combine = TRUE)

datasets <- unique(dc_df[["data"]])
data_dims <- vapply(
  datasets,
  function(x) {
    suppressMessages(
      dim(
        readRDS(paste0("data/", x, ".rds"))
      )
    )
  },
  FUN.VALUE = numeric(2)
)
data_dims <- as.data.frame(t(data_dims))
colnames(data_dims) <- c("nGenes", "nCells")
data_dims[["data"]] <- datasets


source(here("src/analysis/downsampling.R"))
source(here("src/analysis/removing_cells.R"))

source(here("src/analysis/true_positives.R"))

# file_df <- rbind(advi_df, dc_df)
file_df <- dc_df
df <- merge(file_df, data_dims)

source(here("src/analysis/time_plot.R"))

references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)

source(here("src/analysis/de_on_table.R"))
source(here("src/analysis/chain_plots.R"))


source(here("src/analysis/ess.R"))
source(here("src/analysis/hpd.R"))


source(here("src/analysis/batchinfo.R"))
source(here("src/analysis/identifiability.R"))