library("here")

source(here("src/analysis/preamble.R"))
source(here("src/analysis/read_chains.R"))

pes_all <- parallel::mclapply(
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
            mu = colMedians(chain@parameters$mu),
            delta = colMedians(chain@parameters$delta),
            epsilon = colMedians(chain@parameters$epsilon)
        )
    },
    mc.cores = 8
)

pe_all <- bind_rows(pes_all)
pe_all[which(pe_all[["chains"]] == 1), "by"] <- "Reference"
pe_all[pe_all[["by"]] == "advi", "by"] <- "ADVI"
pe_all[pe_all[["by"]] == "gene", "by"] <- "Divide and conquer"

pe_ref <- pe_all[which(pe_all[["by"]] == "Reference"), ]
pe_scal <- pe_all[which(pe_all[["by"]] != "Reference"), ]
pe_sub_by_data <- lapply(unique(pe_scal$data), function(d) {
    pe_sub <- pe_scal[pe_scal$data == d, ]
    pe_sub[, c("mu", "epsilon", "delta")] <- sapply(
        c("mu", "epsilon", "delta"),
        function(p) {
            pe_sub[[p]] - pe_ref[pe_ref$data == d, p]
        }
    )
    pe_sub
})
pe_sub <- do.call(rbind, pe_sub_by_data)

pe_sub[["chains"]] <- factor(
pe_sub[["chains"]],
    levels = sort(unique(pe_sub[["chains"]]))
)

pe_ordered <- pe_sub %>% 
    group_by(data, chains, seed, by) %>% 
    arrange(feature, .by_group = TRUE)

pe_ordered[["data"]] <- sub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    pe_ordered[["data"]],
    perl = TRUE
)

gs <- lapply(c("mu", "delta", "epsilon"),
    function(p) {
        ggplot(pe_ordered) +
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
            # scale_y_log10() +
            labs(
                x = "Number of partitions",
                y = "Difference in point estimates"
            )
    }
)

ggsave(gs[[1]],
    file = "figs/point_estimates_mu.pdf",
    width = 6, height = 4.5
)
ggsave(gs[[2]],
    file = "figs/point_estimates_delta.pdf",
    width = 6, height = 4.5
)
ggsave(gs[[3]],
    file = "figs/point_estimates_epsilon.pdf",
    width = 6, height = 4.5
)
save.image("rdata/point_estimates.RData")
