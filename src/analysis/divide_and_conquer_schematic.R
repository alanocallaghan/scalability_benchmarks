library("ggplot2")
theme_set(theme_classic())

chains_gene <- rbind(
    data.frame(
        chain = "Subposterior 1",
        param = "Gene 1",
        value = rnorm(1000, 3)
    ),
    data.frame(
        chain = "Subposterior 2",
        param = "Gene 2",
        value = rnorm(1000, 5)
    ),
    data.frame(
        chain = "Subposterior 3",
        param = "Gene 3",
        value = rnorm(1000, 8)
    ),
    data.frame(
        chain = "Original posterior",
        param = "Gene 1",
        value = rnorm(1000, 3)
    ),
    data.frame(
        chain = "Original posterior",
        param = "Gene 2",
        value = rnorm(1000, 5)
    ),
    data.frame(
        chain = "Original posterior",
        param = "Gene 3",
        value = rnorm(1000, 8)
    )
)
combined <- chains_gene[1:3000, ]
combined$chain <- "Combined posterior"
chains_gene <- rbind(
    chains_gene, combined
)
chains_cell <- rbind(
    data.frame(
        chain = "Subposterior 1",
        param = "Gene 1",
        value = rnorm(1000, 2.7, 0.7)
    ),
    data.frame(
        chain = "Subposterior 1",
        param = "Gene 2",
        value = rnorm(1000, 5, 1.1)
    ),
    data.frame(
        chain = "Subposterior 1",
        param = "Gene 3",
        value = rnorm(1000, 8)
    ),
    data.frame(
        chain = "Subposterior 2",
        param = "Gene 1",
        value = rnorm(1000, 2.1, 0.9)
    ),
    data.frame(
        chain = "Subposterior 2",
        param = "Gene 2",
        value = rnorm(1000, 5, 0.65)
    ),
    data.frame(
        chain = "Subposterior 2",
        param = "Gene 3",
        value = rnorm(1000, 8, 0.5)
    ),
    data.frame(
        chain = "Subposterior 3",
        param = "Gene 1",
        value = rnorm(1000, 2.5, 0.4)
    ),
    data.frame(
        chain = "Subposterior 3",
        param = "Gene 2",
        value = rnorm(1000, 5, 0.6)
    ),
    data.frame(
        chain = "Subposterior 3",
        param = "Gene 3",
        value = rnorm(1000, 8, 1.4)
    ),
    data.frame(
        chain = "Combined posterior",
        param = "Gene 1",
        value = rnorm(1000, 2.4, 0.8)
    ),
    data.frame(
        chain = "Combined posterior",
        param = "Gene 2",
        value = rnorm(1000, 5, 0.62)
    ),
    data.frame(
        chain = "Combined posterior",
        param = "Gene 3",
        value = rnorm(1000, 8.02, 1.4)
    ),
    data.frame(
        chain = "Original posterior",
        param = "Gene 1",
        value = rnorm(1000, 3, 0.5)
    ),
    data.frame(
        chain = "Original posterior",
        param = "Gene 2",
        value = rnorm(1000, 5)
    ),
    data.frame(
        chain = "Original posterior",
        param = "Gene 3",
        value = rnorm(1000, 8)
    )
)

t <- list(
    aes(x=value, color=param),
    geom_density(),
    facet_wrap(~chain, ncol=1, scale = "free_y"),
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom"),
    scale_color_brewer(palette = "Accent", name = NULL),
    labs(x="Parameter value", y = "Posterior density")
)
levs <- c("Subposterior 1", "Subposterior 2", "Subposterior 3", "Combined posterior", "Original posterior")
chains_gene$chain <- factor(chains_gene$chain, levels = levs)
chains_cell$chain <- factor(chains_cell$chain, levels = levs)

gg <- ggplot(chains_gene) + t + ggtitle("Gene-wise divide and conquer")

gc <- ggplot(chains_cell) + t + ggtitle("Cell-wise divide and conquer")

library("patchwork")
gg + gc + plot_layout(guides = "collect") +
    plot_annotation(tag_level="A") &
    theme(legend.position='bottom')

ggsave("figs/divide_and_conquer_schematic.pdf", width = 7, height = 6)
