library("SingleCellExperiment")
library("ggplot2")
library("cowplot")

source("src/analysis/functions.R")
theme_set(theme_bw())

datasets <- c(
    "tung",
    "buettner",
    # "pbmc",
    "ibarra-soria",
    "chen",
    "zeisel"
)
datas <- lapply(datasets, function(d) readRDS(paste0("rdata/", d, ".rds")))

names(datas) <- c(
    "Tung",
    "Buettner",
    # "10x PBMC",
    "Ibarra-Soria",
    "Chen",
    "Zeisel"
)

l <- lapply(names(datas),
    function(n) {
        x <- datas[[n]]
        data.frame(
            name = n,
            lib_sizes = colSums(counts(x)),
            num_feats = colSums(counts(x) != 0)
        )
    }
)
data_df <- do.call(rbind, l)


scale <- scale_color_brewer(
    palette = "Set2",
    name = "Dataset    ",
    aesthetics = c("color", "fill")
)


g1 <- ggplot(data_df, aes(x = lib_sizes, color = name, fill = name)) +
    labs(x = "Total counts per cell", y = "Density") +
    scale +
    geom_density(alpha = 0.2) +
    theme(
        legend.position = "bottom",
        panel.grid = element_blank()
    ) +
    scale_x_log10()
ggsave("figs/libsize_density.pdf", width = 5, height = 4)


g2 <- ggplot(data_df, aes(x = num_feats, color = name, fill = name)) +
    labs(x = "Number of expressed features per cell", y = "Density") +
    scale +
    geom_density(alpha = 0.2) +
    theme(
        legend.position = "bottom",
        panel.grid = element_blank()
    ) +
    scale_x_log10()
ggsave("figs/complexity_density.pdf", width = 5, height = 4)

ggp <- plot_with_legend_below(g1, g2, rel_heights = c(0.9, 0.1))
ggsave("figs/cell_plots.pdf", width = 7, height = 4)


l <- lapply(names(datas),
    function(n) {
      x <- datas[[n]]
      data.frame(
          name = n,
          mean_expression = rowMeans(counts(x)),
          dropout = rowMeans(counts(x) == 0)
      )
    }
)
data_df <- do.call(rbind, l)

g3 <- ggplot(data_df, aes(x = mean_expression, color = name, fill = name)) +
    labs(x = "Mean count per gene", y = "Density") +
    scale +
    geom_density(alpha = 0.2) +
    theme(
        legend.position = "bottom",
        panel.grid = element_blank()
    ) +
    scale_x_log10()
ggsave("figs/expression_density.pdf", width = 5, height = 4)

g4 <- ggplot(data_df, aes(x = dropout, color = name, fill = name)) +
    labs(x = "Proportion of zeros per gene", y = "Density") +
    scale +
    theme(
        legend.position = "bottom",
        panel.grid = element_blank()
    ) +
    geom_density(alpha = 0.2)

ggsave("figs/dropout_density.pdf", width = 5, height = 4)

ggp <- plot_with_legend_below(g1, g2, rel_heights = c(0.9, 0.1))
ggsave("figs/gene_plots.pdf", width = 7, height = 4)


library("patchwork")
(g1 + g2) / (g3 + g4) + plot_layout(guides="collect") +
    plot_annotation(tag_level="A") &
    theme(legend.position="bottom")

ggsave("figs/data_plots.pdf", width = 7, height = 5)
