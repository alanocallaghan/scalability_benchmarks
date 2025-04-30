suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("ggbeeswarm")
  library("here")
  library("BASiCS")
  library("coda")
  library("viridis")
})

source(here("src/analysis/functions.R"))
theme_set(theme_bw())


dir.create("figs/elbo", recursive = TRUE, showWarnings = FALSE)

advi_files <- list.files("outputs/advi", full.names = TRUE)
advi_triplets <- file2triplets(advi_files)
advi_elbo <- lapply(advi_triplets, function(x) readRDS(x[[3]]))
advi_triplets <- lapply(advi_triplets, function(x) x[-3])
advi_df <- read_triplets(advi_triplets)

parsed_elbos <- lapply(advi_elbo, parse_elbo)

elbo_df <- as.data.frame(advi_df)
elbo_df$elbos <- parsed_elbos
elbo_df$elbos <- lapply(seq_len(nrow(elbo_df)), 
    function(i) {
        d <- elbo_df$elbos[[i]]
        d$data <- elbo_df[[i, "data"]]
        d$seed <- elbo_df[[i, "seed"]]
        d
    }
)

elbo_df <- do.call(rbind, elbo_df$elbos)
elbo_df$Data <- sub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    elbo_df$data,
    perl = TRUE
)

plots <- lapply(unique(elbo_df$data),
    function(d) {
        df <- elbo_df[elbo_df$data == d, ]
        D <- df$Data[[1]]
        g <- ggplot(df) +
            aes(x = iter, y = ELBO, color = factor(seed)) +
            geom_line() +
            theme(legend.position = "none") +
            # scale_y_log10() +
            # ggtitle(D) +
            labs(x = "Iteration") +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin = unit(c(0, 0.05, 0, 0.025), "npc")
            )
            labs(x = "Iteration")
        ggsave(g, file = sprintf("figs/elbo_%s.pdf", d), width = 4, height = 3)
        g
    }
)
