###############################################################################
##
## Run TestDE
##
###############################################################################

df <- do_de(df, ref_df = references, match_column = "data", data_dims,
    mc.cores = 12)


###############################################################################
##
## DE numbers
##
###############################################################################

sdf_de  <- df[,
    c("data", "chains", "pDiffExp", "pDiffDisp", "pDiffResDisp")
]
mdf_de <- reshape2::melt(sdf_de, id.vars = c("data", "chains"))
mdf_de$chains <- as.numeric(mdf_de$chains)
mdf_de$variable <- gsub("^pDiffExp$", "mu", mdf_de$variable)
mdf_de$variable <- gsub("^pDiffDisp$", "delta", mdf_de$variable)
mdf_de$variable <- gsub("^pDiffResDisp$", "epsilon", mdf_de$variable)
mdf_de$variable <- factor(
    mdf_de$variable,
    levels = c("mu", "delta", "epsilon")
)


mdf_de$data <- sub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    mdf_de$data,
    perl = TRUE
)
# mdf_de$data_desc <- paste(mdf_de$data, "data:\n",
#     data_dims[tolower(mdf_de$data), "nGenes"], "genes,",
#     data_dims[tolower(mdf_de$data), "nCells"], "cells"
# )
mdf_de$data_desc <- paste(mdf_de$data, "et al. data")
mdf_de <- mdf_de[mdf_de$variable %in% c("mu", "epsilon"), ]
mdf_de$variable <- droplevels(mdf_de$variable)
levels(mdf_de$variable) <- c("Mean", "Residual over-dispersion")


# mdf_de <- mdf_de %>% filter(variable %in% c("mu", "epsilon"))

advi_mdf_de <- mdf_de %>%
    filter(is.na(chains)) %>%
    group_by(data_desc, variable, chains) %>%
    summarize(value = median(value))


g <- ggplot(mdf_de[!(is.na(mdf_de$chains) | mdf_de$chains == 1), ],
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
        name = "Spuriously differentially expressed genes",
        labels = scales::percent
    ) +
    theme(
        # text = element_text(size = 18),
        legend.position = "bottom",
        panel.grid = element_blank()
    ) +
    scale_color_brewer(name = "Parameter", palette = "Dark2")

ggsave(here("figs/diffexp_plot.pdf"), width = 5, height = 5)
save.image("rdata/diffexp_plot.RData")
