library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")

source(here("src/analysis/preamble.R"))

theme_set(theme_bw())
source(here("src/analysis/functions.R"))

ds_files <- list.files("outputs/downsampling/divide/", full.names = TRUE)
dst <- file2triplets(ds_files)
dst <- dst[as.logical(sapply(dst, length))]
ds_df <- read_triplets(dst, combine = TRUE)

ref_files <- list.files("outputs/downsampling/reference", full.names = TRUE)
ref_df_ds <- read_triplets(file2triplets(ref_files), combine = TRUE)
ref_df_ds$chain <- lapply(ref_df_ds$file, readRDS)
ds_df <- do_de(ds_df, ref_df_ds, "downsample_rate", data_dims, mc.cores = 12)

mdf_ds <- reshape2::melt(ds_df,
    measure.vars = c("pDiffExp", "pDiffDisp", "pDiffResDisp")
)
mdf_ds$variable <- gsub("pDiffExp", "mu", mdf_ds$variable)
mdf_ds$variable <- gsub("pDiffDisp", "delta", mdf_ds$variable)
mdf_ds$variable <- gsub("pDiffResDisp", "epsilon", mdf_ds$variable)
mdf_ds$variable <- factor(mdf_ds$variable, levels = c("mu", "delta", "epsilon"))
mdf_ds <- merge(mdf_ds, data_dims)


# ref_df_ds[[1, "data"]]
# sce <- readRDS(paste0("rdata/", "tung", ".rds"))
# libsize <- colSums(counts(sce))
mdf_ds$mean_libsize <- mdf_ds$libsize * mdf_ds$downsample_rate
mdf_ds$downsample_rate_t <- factor(
    paste(mdf_ds$downsample_rate * 100, "%"),
    levels = paste0(
        sort(unique(ds_df$downsample_rate), decreasing = TRUE) * 100,
        "%"
    )
)

mdf_ds_sub <- mdf_ds[mdf_ds$data == "chen", ]
mdf_ds_sub <- mdf_ds_sub[mdf_ds_sub$variable %in% c("mu", "epsilon"), ]
mdf_ds_sub$variable <- droplevels(mdf_ds_sub$variable)
levels(mdf_ds_sub$variable) <- c("Mean", "Residual over-dispersion")


g <- ggplot(mdf_ds_sub) +
    aes(
        x = factor(format(signif(mean_libsize, digits = 3), big.mark = ",")),
        y = value,
        color = variable
    ) +
    geom_quasirandom(dodge.width = 0.25, size = 0.7, groupOnX = TRUE) +
    # facet_wrap(~data) +
    scale_color_brewer(name = "Parameter", palette = "Dark2") +
    scale_y_continuous(
        label = scales::percent,
        limits = c(0, max(0.25, max(mdf_ds_sub$value)))
    ) +
    theme(legend.position = "bottom", panel.grid = element_blank()) +
    labs(
        x = "Expected median library size",
        y = "Spuriously differentially expressed genes"
    )

ggsave("figs/downsampling.pdf", width = 5, height = 4)
save.image("rdata/downsampling.RData")
