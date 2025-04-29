library("BASiCS")
library("ggplot2")
library("here")
library("reshape2")

theme_set(theme_bw())


source(here("src/analysis/preamble.R"))
dc_files <- list.files("outputs/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(file2triplets(dc_files), combine = TRUE)

file_df <- dc_df
df <- merge(file_df, data_dims)

references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)

cs_files <- list.files("outputs/cell_splitting", full.names = TRUE)
cs_df <- data.frame(
    file = cs_files,
    data = gsub(".*/(\\w+)_chains-\\d+_seed-\\d+\\.rds", "\\1", cs_files),
    chains = gsub(".*/\\w+_chains-(\\d+)_seed-\\d+\\.rds", "\\1", cs_files),
    by = "cell"
)
cs_df <- merge(cs_df, data_dims)

cs_de_df <- do_de(cs_df, ref_df = references, match_column = "data", data_dims,
  mc.cores = 2)

cs_de_df$chain <- NULL
cs_de_df$data <- gsub("([\\w])([\\w]+)", "\\U\\1\\L\\2", cs_de_df$data, perl = TRUE)
cs_de_df_sub <- cs_de_df[, c("data", "chains", "pDiffExp", "pDiffDisp", "pDiffResDisp")]
mdf_cs <- melt(cs_de_df_sub, id.vars = c("data", "chains"))
mdf_cs$chains <- factor(mdf_cs$chains, levels = c(2, 4, 8, 16))
levels(mdf_cs$variable) <- c("mu", "delta", "epsilon")


g <- ggplot(mdf_cs) +
    aes(x = chains, colour = variable, y = value) +
    geom_quasirandom(groupOnX = TRUE, dodge.width = 0.5, size = 0.5) +
    facet_wrap(~data) +
    scale_colour_brewer(
        name = "Parameter",
        palette = "Set1"
    ) +
    labs(x = "Number of partitions") +
    scale_y_continuous(
        name = "Proportion of genes differentially expressed",
        labels = scales::percent
    ) +
    theme(
        legend.position = "bottom",
        panel.grid = element_blank()
    )

ggsave("figs/cell_splitting.pdf", width = 5, height = 5)
save.image("rdata/cell_splitting.RData")
