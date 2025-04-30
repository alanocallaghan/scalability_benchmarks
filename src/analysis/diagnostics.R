library("here")
library("BASiCS")
source(here("src/analysis/preamble.R"))
source(here("src/analysis/read_chains.R"))

ess <- function (x) {
    vars <- matrixStats::colVars(x)
    spec <- numeric(ncol(x))
    has_var <- vars > 1e-10
    if (any(has_var, na.rm = TRUE)) {
        spec[which(has_var)] <- apply(x[, which(has_var), drop = FALSE],
            2, function(y) {
                a <- ar(y, aic = TRUE)
                a$var.pred / (1 - sum(a$ar))^2
            })
    }
    setNames(ifelse(spec == 0, 0, nrow(x) * vars/spec), colnames(x))
}

for (measure in c("ess", "geweke.diag")) {

    dir.create(
        sprintf("figs/%s", gsub("\\.", "_", measure)),
        showWarnings = FALSE
    )

    measure_name <- BASiCS:::.MeasureName(measure)
    diag_all_list <- parallel::mclapply(
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
                mu = BASiCS:::.GetMeasure(chain, "mu", measure),
                delta = BASiCS:::.GetMeasure(chain, "delta", measure),
                epsilon = BASiCS:::.GetMeasure(chain, "epsilon", measure)
            )
        }, mc.cores = getOption("mc.cores", 4)
    )
    diag_all <- bind_rows(diag_all_list)
    diag_all[which(diag_all[["chains"]] == 1), "by"] <- "Reference"
    diag_all <- diag_all[diag_all[["by"]] != "advi", ]
    diag_all[diag_all[["by"]] == "gene", "by"] <- "Divide and conquer"
    diag_all[["chains"]] <- factor(
        diag_all[["chains"]],
        levels = sort(unique(diag_all[["chains"]]))
    )

    diag_all[["data"]] <- sub(
        "([\\w])([\\w]+)", "\\U\\1\\L\\2",
        diag_all[["data"]],
        perl = TRUE
    )

    scale <- if (measure == "ess") {
        list(
            scale_y_log10(name = measure_name),
            geom_hline(
                yintercept = 100,
                col = "firebrick",
                linetype = "dashed"
            )
        )
    } else {
        list(
            scale_y_continuous(name = measure_name),
            geom_hline(
                yintercept = c(-3, 3),
                col = "firebrick",
                linetype = "dashed"
            )
        )
    }

    gs <- lapply(c("mu", "delta", "epsilon"),
        function(param) {
            ggplot(diag_all) +
                aes_string(x = "chains", y = param, color = "by", fill = "by") +
                # geom_jitter(height = 0, width = 0.2) +
                geom_violin(alpha = 0.2) +
                geom_boxplot(alpha = 0.2, width = 0.1, outlier.colour = NA) +
                facet_wrap(~data, scales = "free_y") +
                scale_fill_brewer(
                    name = "Inference method",
                    palette = "Dark2",
                    aesthetics = c("fill", "color")
                ) +
                scale +
                theme(
                    panel.grid = element_blank(),
                    legend.position = "bottom"
                ) +
                labs(
                    x = "Number of partitions",
                    y = "Effective sample size"
                )
        }
    )
    names(gs) <- c("mu", "delta", "epsilon")
    for (param in names(gs)) {
        ggsave(gs[[param]],
            file = sprintf("figs/%s_%s_all.pdf", gsub("\\.", "_", measure), param),
            width = 5, height = 4
        )
    }
}
