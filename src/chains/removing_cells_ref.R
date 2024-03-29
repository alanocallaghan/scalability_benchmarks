#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-f", "--fraction", type = "double")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")

args <- parser$parse_args()

source(here("src/chains/benchmark_code.R"))

set.seed(42)
data <- readRDS(here("rdata", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

frac <- args[["fraction"]]
data <- data[, sample(ncol(counts(data)), floor(ncol(counts(data)) * frac))]
if (length(altExpNames(data))) {
    spikes <- altExp(data, "spike-ins")
    ind_keep_spike <- rowSums(assay(spikes)) != 0
    metadata(data)$SpikeInput <- metadata(data)$SpikeInput[ind_keep_spike, ]
    altExp(data, "spike-ins") <- spikes[ind_keep_spike, ]
}


time <- system.time({
    chain <- BASiCS_MCMC(
        data,
        Regression = TRUE,
        WithSpikes = as.logical(length(altExpNames(data))),
        PrintProgress = FALSE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )
})
config <- list(
    data = args[["data"]],
    chains = 1,
    by = "gene",
    seed = 42,
    proportion_retained = args[["fraction"]]
)
saveRDS(chain, file = file.path(dir, "chains.rds"))
saveRDS(time, file = file.path(dir, "time.rds"))
saveRDS(config, file = file.path(dir, "config.rds"))
