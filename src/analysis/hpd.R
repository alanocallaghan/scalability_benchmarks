library("here")
source(here("src/analysis/preamble.R"))
source(here("src/analysis/read_chains.R"))

get_hpd_width <- function(chain, param) {
    hpd <- HPDinterval(as.mcmc(chain@parameters[[param]]))
    (hpd[, "upper"] - hpd[, "lower"])
}

hpds_all <- parallel::mclapply(
    seq_len(nrow(df)),
    function(i) {
        cat(i, "/", nrow(df), "\n")
        chain <- readRDS(df[[i, "file"]])
        data.frame(
            data = df[[i, "data"]],
            chains = df[[i, "chains"]],
            seed = df[[i, "seed"]],
            by = df[[i, "by"]],
            feature = colnames(chain@parameters[["mu"]]),
            mu = get_hpd_width(chain, "mu"),
            delta = get_hpd_width(chain, "delta"),
            epsilon = get_hpd_width(chain, "epsilon")
        )
    }, mc.cores = getOption("mc.cores", 6)
)

hpdf_all <- bind_rows(hpds_all)
hpdf_all[which(hpdf_all[["chains"]] == 1), "by"] <- "Reference"
hpdf_all[hpdf_all[["by"]] == "advi", "by"] <- "ADVI"
hpdf_all[hpdf_all[["by"]] == "gene", "by"] <- "Divide and conquer"
hpdf_all[["chains"]] <- factor(
    hpdf_all[["chains"]],
    levels = sort(unique(hpdf_all[["chains"]]))
)

hpdf_ordered <- hpdf_all %>% 
    group_by(data, chains, seed, by) %>% 
    arrange(feature, .by_group = TRUE)

hpdf_ordered[["data"]] <- sub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    hpdf_ordered[["data"]],
    perl = TRUE
)

gs <- lapply(c("mu", "delta", "epsilon"),
    function(p) {
        ggplot(hpdf_ordered) +
            aes_string(y = p, x = "chains", color = "by", fill = "by") +
            geom_violin(alpha = 0.2) +
            geom_boxplot(alpha = 0.2, width = 0.1, outlier.colour = NA) +
            facet_wrap(~data, scales = "free_y") +
            scale_fill_brewer(
                name = "Inference method",
                palette = "Dark2",
                aesthetics = c("fill", "color")
            ) +
            theme(
                legend.position = "bottom",
                panel.grid = element_blank()
            ) +
            scale_y_log10() +
            labs(
                x = "Number of partitions",
                y = "95% HPD interval width"
            )
    }
)

ggsave(gs[[1]],
    file = "figs/hpd_width_mu.pdf",
    width = 6, height = 5
)
ggsave(gs[[2]],
    file = "figs/hpd_width_delta.pdf",
    width = 6, height = 5
)
ggsave(gs[[3]],
    file = "figs/hpd_width_epsilon.pdf",
    width = 6, height = 5
)
