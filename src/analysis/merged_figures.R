# merging figures
library("ggplot2") ## 3.4.0!!!!
library("patchwork")
library("ggbeeswarm")
library("parallel")
theme_set(theme_bw())
files <- c(
    "diffexp_plot.RData",
    "true_positives.RData",
    "downsampling.RData",
    "removing_cells.RData"
)

envs <- mclapply(files, function(file) {
    env <- new.env()
    load(file.path("rdata", file), envir = env)
    env
}, mc.cores=4)

gs <- lapply(envs, function(env) {
    env$g
})

advi_mdf_de <- envs[[1]]$advi_mdf_de
mdf_de <- envs[[1]]$mdf_de
gs[[1]] <- ggplot(mdf_de[!(is.na(mdf_de$chains) | mdf_de$chains == 1), ],
      aes(
          x = factor(chains),
          y = value,
          # group = data,
          color = variable
      )
    ) +
    geom_hline(
        data = advi_mdf_de,
        alpha = 0.5,
        linetype = "dashed",
        show_guide=FALSE,
        aes(
            yintercept = value,
            color = variable
        )
    ) +
    geom_quasirandom(
        groupOnX = TRUE,
        dodge.width = 0.5,
        size = 0.7
    ) +
    facet_wrap(~data_desc, nrow = 2, ncol = 2, scales = "free_y") +
    scale_x_discrete(name = "Number of subsets") +
    scale_y_continuous(
        name = "Portion of genes differentially expressed",
        labels = scales::percent
    ) +
    theme(
        # text = element_text(size = 18),
        legend.position = "bottom",
        panel.grid = element_blank()
    ) +
    scale_color_brewer(name = "Parameter", palette = "Dark2")

gs <- lapply(gs, function(g) g + scale_color_brewer(name = "Parameter", palette = "Dark2"))

wrap_elements(plot=gs[[1]] + guides(colour = "none")) + gs[[2]] + gs[[3]] + gs[[4]] +
    plot_layout(guides="collect") +
    plot_annotation(tag_level="A") &
    theme(legend.position="bottom")



ggsave("figs/accuracy_plots.pdf", width=9, height=9)
