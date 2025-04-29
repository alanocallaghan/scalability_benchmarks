library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")

source(here("src/analysis/preamble.R"))

theme_set(theme_bw())
source(here("src/analysis/functions.R"))


rm_files <- list.files("outputs/removing/divide", full.names = TRUE)
rm_files <- rm_files[grep("chen", rm_files)]
rm_ref_files <- list.files("outputs/removing/reference", full.names = TRUE)
rm_ref_files <- rm_ref_files[grep("chen", rm_ref_files)]


rmt <- file2triplets(rm_files)
rmt <- rmt[as.logical(sapply(rmt, length))]
rm_df <- read_triplets(rmt, combine = TRUE)

ref_df_rm <- read_triplets(file2triplets(rm_ref_files), combine = TRUE)
ref_df_rm$chain <- lapply(ref_df_rm$file, readRDS)

rm_df <- do_de(
    rm_df,
    ref_df = ref_df_rm,
    match_column = "proportion_retained",
    data_dims,
    mc.cores = 12
)

mdf_rm <- reshape2::melt(
    rm_df,
    measure.vars = c("pDiffExp", "pDiffDisp", "pDiffResDisp")
)
mdf_rm$variable <- gsub("pDiffExp", "mu", mdf_rm$variable)
mdf_rm$variable <- gsub("pDiffDisp", "delta", mdf_rm$variable)
mdf_rm$variable <- gsub("pDiffResDisp", "epsilon", mdf_rm$variable)
mdf_rm$variable <- factor(mdf_rm$variable, levels = c("mu", "delta", "epsilon"))
mdf_rm$cells_retained <- mdf_rm$proportion_retained * mdf_rm$nCells

mdf_rm$proportion_retained <- factor(
    paste(mdf_rm$proportion_retained * 100, "%"),
    levels = paste(
        sort(unique(rm_df$proportion_retained), decreasing = TRUE) * 100,
        "%"
    )
)

# mdf_rm_sub <- mdf_rm[mdf_rm$data == "zeisel", ]
mdf_rm_sub <- mdf_rm[mdf_rm$data == "chen", ]
mdf_rm_sub <- mdf_rm_sub[mdf_rm_sub$variable %in% c("mu", "epsilon"), ]
mdf_rm_sub$variable <- droplevels(mdf_rm_sub$variable)
levels(mdf_rm_sub$variable) <- c("Mean", "Residual over-dispersion")


g <- ggplot(mdf_rm_sub) +
    aes(x = factor(round(cells_retained)), y = value, color = variable) +
    geom_quasirandom(dodge.width = 0.25, size = 0.7, groupOnX = TRUE) +
    scale_color_brewer(name = "Parameter", palette = "Dark2") +
    # scale_x_reverse("Number of cells") +
    # facet_wrap(~data) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        x = "Number of cells in dataset",
        y = "Spuriously differentially expressed genes"
    ) +
    theme(
        axis.text.x = element_text(hjust = 1, angle = 45),
        panel.grid = element_blank(),
        legend.position = "bottom"
    )

ggsave("figs/removing_cells.pdf", width = 5, height = 4)
save.image("rdata/removing_cells.RData")
