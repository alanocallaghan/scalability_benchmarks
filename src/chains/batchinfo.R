#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


data <- readRDS(paste0("data/", args[["data"]], ".rds"))

N <- 20000
Thin <- 10
Burn <- 10000

bi <- BASiCS_MCMC(
  data,
  N = args[["iterations"]],
  Thin = max((args[["iterations"]] / 2) / 1000, 2),
  Burn = max(args[["iterations"]] / 2, 4),
  PrintProgress = FALSE,
  WithSpikes = FALSE,
  Regression = TRUE
)

data@colData$BatchInfo <- 1

non_bi <- BASiCS_MCMC(
  data,
  N = N,
  Thin = Thin,
  Burn = Burn,
  PrintProgress = FALSE,
  WithSpikes = FALSE,
  Regression = TRUE
)

dir <- args[["output"]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(bi, file.path(dir, "batch.rds"))
saveRDS(non_bi, file.path(dir, "/nobatch.rds"))


library(ggplot2)
theme_set(theme_bw())
x <- rgamma(10000, 1, 1)
ggplot() +
  aes(x) +
  geom_histogram(boundary=0, bins=nclass.FD(x)) +
  geom_vline(aes(xintercept = quantile(x)), col="blue", alpha = 0.7) +
  geom_vline(aes(xintercept = seq(min(x), max(x), length.out=5)), col="red", alpha = 0.7)