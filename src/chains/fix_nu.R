#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-d", "--dataset", type = "character")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

dataset <- args[["dataset"]]
sce <- readRDS(paste0("rdata/", dataset, ".rds"))

if (dataset == "chen") {
    fit_fix <- BASiCS_MCMC(
        sce,
        PrintProgress = FALSE,
        FixNu = TRUE,
        WithSpikes = FALSE,
        Regression = TRUE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )

    fit_var <- BASiCS_MCMC(
        sce,
        PrintProgress = FALSE,
        FixNu = FALSE,
        WithSpikes = FALSE,
        Regression = TRUE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )
} else {

    droplet_sce <- readRDS("rdata/ibarra-soria.rds")
    ind_presom <- colData(droplet_sce)[["Cell_type"]] == "PSM"
    ind_som <- colData(droplet_sce)[["Cell_type"]] == "SM"


    PSM_Data <- droplet_sce[, ind_presom]
    SM_Data <- droplet_sce[, ind_som]
    
    fit_fix <- list()
    fit_var <- list()

    fit_fix[["PSM"]] <- BASiCS_MCMC(
        PSM_Data,
        PrintProgress = FALSE,
        FixNu = TRUE,
        WithSpikes = FALSE,
        Regression = TRUE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )

    fit_var[["PSM"]] <- BASiCS_MCMC(
        PSM_Data,
        PrintProgress = FALSE,
        FixNu = FALSE,
        WithSpikes = FALSE,
        Regression = TRUE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )
    fit_fix[["SM"]] <- BASiCS_MCMC(
        SM_Data,
        PrintProgress = FALSE,
        FixNu = TRUE,
        WithSpikes = FALSE,
        Regression = TRUE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )

    fit_var[["SM"]] <- BASiCS_MCMC(
        SM_Data,
        PrintProgress = FALSE,
        FixNu = FALSE,
        WithSpikes = FALSE,
        Regression = TRUE,
        N = args[["iterations"]],
        Thin = max((args[["iterations"]] / 2) / 1000, 2),
        Burn = max(args[["iterations"]] / 2, 4)
    )

}

dir.create("outputs/fix_nu", showWarnings = FALSE)
saveRDS(fit_fix, file.path(args[["output"]], paste0(dataset, "-fix.rds")))
saveRDS(fit_var, file.path(args[["output"]], paste0(dataset, "-var.rds")))
