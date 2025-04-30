read_triplets <- function(triplets, combine = FALSE) {
    rows <- lapply(
        triplets,
        function(x) {
            row <- readRDS(x[[2]])
            row <- as.data.frame(row)
            row[["time"]] <- readRDS(x[[3]])[["elapsed"]]
            row
        }
    )
    df <- do.call(rbind, rows)
    df[["file"]] <- lapply(triplets, function(x) x[[1]])
    df
}

file2triplets <- function(files) {
  lapply(files, list.files, full.names = TRUE)
}

## Convergence first
parse_elbo <- function(x) {
    x <- gsub("Chain 1:\\s+", "", x)
    normal <- grep("Drawing", x)
    abnormal <- grep("Informational", x)
    if (!length(normal)) {
        stop("Something is deeply wrong here")
    }
    ind_end <- normal - 2
    if (length(abnormal)) {
        ind_end <- abnormal - 1
    }
    elbo <- x[(grep("Begin stochastic", x) + 1):ind_end]
    ## This ain't quite right
    elbo <- gsub("MAY BE DIVERGING... INSPECT ELBO", "", elbo, fixed = TRUE)
    elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
        elbo[-c(1, grep("CONVERGED", elbo))],
        "NOTCONVERGED"
    )
    ## If both mean and median elbo converge this ends up with CONVERGED CONVERGED
    ## and therefore another column
    elbo <- gsub("(MEDIAN |MEAN )?ELBO CONVERGED", "CONVERGED", elbo)
    elbo <- gsub("(\\s+CONVERGED){2}", " CONVERGED", elbo)
    elbo <- strsplit(elbo, "\\s+")
    elbo <- do.call(rbind, elbo)
    colnames(elbo) <- elbo[1, ]
    elbo <- elbo[-1, ]
    elbo <- as.data.frame(elbo, stringsAsFactors = FALSE)
    elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
    elbo
}

do_de <- function(
        df,
        ref_df,
        match_column,
        data_dims,
        mc.cores = getOption("mc.cores", 2)
    ) {

    edr <- parallel::mclapply(
        seq_len(nrow(df)),
        function(i) {
        cat(i, "/", nrow(df), "\n")
        if (isTRUE(df[[i, "chains"]] == 1)) {
            return(rep(list(NULL), 3))
        }
        ind <- ref_df[[match_column]] == df[[i, match_column]] &
            ref_df[["data"]] == df[i, "data"]
        chain <- readRDS(df[[i, "file"]])
        if (length(chain) > 1) {
            suppressMessages(
            chain <- BASiCS:::.combine_subposteriors(
                chain,
                SubsetBy = "gene"
            )
            )
        }
        cp <- intersect(c("nu", "s", "phi"), names(chain@parameters))
        chain@parameters[cp] <- lapply(
            chain@parameters[cp],
            function(x) {
            colnames(x) <- gsub("_Batch.*", "", colnames(x))
            x
            }
        )
        chain@parameters <- BASiCS:::.reorder_params(
            chain@parameters,
            GeneOrder = rownames(ref_df[[which(ind)[[1]], "chain"]])
        )
        nsamples <- nrow(
            ref_df[[which(ind)[[1]], "chain"]]@parameters[["mu"]]
        )
        params <- setdiff(names(chain@parameters), "RBFLocations")
        chain@parameters[params] <- lapply(
            chain@parameters[params],
            function(x) {
            x[seq_len(nsamples), ]
            }
        )
        suppressMessages(
            de <- BASiCS_TestDE(
            ref_df[[which(ind)[[1]], "chain"]],
            chain,
            Plot = FALSE,
            PlotOffset = FALSE,
            EFDR_M = NULL,
            EFDR_D = NULL,
            EFDR_R = NULL,
            MinESS = NA
            )
        )
        lapply(
            de@Results,
            function(x) {
            l <- as.data.frame(x)[["GeneName"]]
            if (!length(l)) NULL else l
            }
        )
        }, mc.cores = mc.cores
    )

    edr_df <- do.call(rbind, edr)
    colnames(edr_df) <- c("DiffExp", "DiffDisp", "DiffResDisp")
    edr_df <- as.data.frame(edr_df)
    df <- cbind(df, edr_df)

    df[, c("nDiffExp", "nDiffDisp", "nDiffResDisp")] <- lapply(
        df[, c("DiffExp", "DiffDisp", "DiffResDisp")],
        function(x) sapply(x, length)
    )

    df <- merge(df, data_dims)
    df[c("pDiffExp", "pDiffDisp", "pDiffResDisp")] <- round(
        df[c("nDiffExp", "nDiffDisp", "nDiffResDisp")] / df[["nGenes"]],
        digits = 3
    )
    df
}

jaccard <- function(a, b) {
    a <- na.omit(a)
    b <- na.omit(b)
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    if (union == 0) 0 else intersection / union
}

sourceme <- function(x) {
    cat("Running", x, "\n")
    source(x)
}

do_fit_plot <- function(i, maxdf, references) {
    cat(i, "/", nrow(maxdf), "\n")
    maxdf <- maxdf[i, ]
    cc <- readRDS(maxdf[["file"]][[1]])
    b <- maxdf[["by"]]
    if (is.list(cc)) {  
        suppressMessages(
        cc <- BASiCS:::.combine_subposteriors(
            cc,
            SubsetBy = "gene"
        )
        )
    }
    d <- maxdf[["data"]]
    rc <- references[[which(references$data == d)[[1]], "chain"]]
    l2 <- if (b == "advi") "ADVI" else "Divide and conquer"
    cc <- BASiCS:::.offset_correct(cc, rc)
    g1 <- BASiCS_ShowFit(rc)
    g2 <- BASiCS_ShowFit(cc)
    ggsave(g1,
        file = here(sprintf("figs/fit/%s_ref.pdf", d)),
        width = 4,
        height = 3
    )
    ggsave(g2,
        file = here(sprintf("figs/fit/%s_%s.pdf", d, b)),
        width = 4,
        height = 3
    )
    cowplot::plot_grid(
        g1 + ggtitle("Reference"),
        g2 + ggtitle(l2)
    )
}

do_de_plot <- function(i, maxdf, references) {
    cat(i, "/", nrow(maxdf), "\n")
    maxdf <- maxdf[i, ]
    d <- maxdf[["data"]]
    b <- maxdf[["by"]]
    rc <- references[[which(references$data == d)[[1]], "chain"]]
    cc <- readRDS(maxdf[["file"]][[1]])
    if (is.list(cc)) {
        suppressMessages(
        cc <- BASiCS:::.combine_subposteriors(
            cc,
            GeneOrder = rownames(rc),
            CellOrder = colnames(rc),
            SubsetBy = "gene"
        )
        )
    }
    cc@parameters <- BASiCS:::.reorder_params(
        cc@parameters,
        GeneOrder = rownames(rc)
    )

    l2 <- if (b == "advi") "ADVI" else "divide & conquer"
    suppressMessages(
        de <- BASiCS_TestDE(
        rc, cc,
        GroupLabel1 = "reference",
        GroupLabel2 = l2,
        Plot = FALSE,
        PlotOffset = FALSE,
        EFDR_M = NULL,
        EFDR_D = NULL,
        EFDR_R = NULL,
        MinESS = NA
        )
    )
    g1 <- BASiCS_PlotDE(de@Results[[1]],
        Plots = "MA",
        Mu = de@Results$Mean@Table$MeanOverall
    ) +
        theme(legend.position = "bottom") +
        ylab("log2(fold change)\nreference vs divide and conquer")
    g2 <- BASiCS_PlotDE(de@Results[[1]],
        Plots = "Volcano",
        Mu = de@Results$Mean@Table$MeanOverall
    ) +
        xlab("log2(fold change)\nreference vs divide and conquer") +
        theme(legend.position = "bottom")
    all <- plot_with_legend_below(g1, g2)
    ggsave(
        file = here(sprintf("figs/de/mu_%s_%s.pdf", d, b)),
        width = 8, height = 3
    )
    g1 <- BASiCS_PlotDE(de@Results[[3]],
        Plots = "MA",
        Mu = de@Results$Mean@Table$MeanOverall
    ) +
        theme(legend.position = "bottom") +
        ylab("Difference\nreference vs divide and conquer")
    g2 <- BASiCS_PlotDE(de@Results[[3]],
        Plots = "Volcano",
        Mu = de@Results$Mean@Table$MeanOverall
    ) +
        theme(legend.position = "bottom") +
        xlab("Difference\nreference vs divide and conquer")
    all <- plot_with_legend_below(g1, g2)
    ggsave(all,
        file = here(sprintf("figs/de/epsilon_%s_%s.pdf", d, b)),
        width = 8, height = 3
    )
    g
}

plot_hpd_interval_width <- function(
        chain1,
        xname,
        chain2,
        yname,
        param,
        type = c("ma", "std"),
        log = FALSE,
        ...
    ) {

    x <- extract_hpd_interval_width(chain1, param)
    y <- extract_hpd_interval_width(chain2, param)
    df <- data.frame(
        x,
        y
    )

    type <- match.arg(type)
    if (type == "ma") {
        if (log) {
        df[] <- lapply(df, log10)
        }
        fun <- param_plot_ma
        mu_df <- data.frame(
        colMeans(chain1@parameters[["mu"]]),
        colMeans(chain2@parameters[["mu"]])
        )
        mus <- rowMeans(mu_df)
    } else {
        fun <- param_plot_std
        mus <- NULL
    }

    colnames(df) <- c(xname, yname)
    g <- fun(
        df,
        xname = xname,
        yname = yname,
        title = paste(
        if (type == "ma") "MA plot" else "Scatter plot",
        "of 95% HPD interval\nof", param, "parameter in",
        yname, "vs", xname, "MCMC"
        ),
        log = log && type == "std",
        mus = mus,
        ...
    )

    # ggMarginal(g)
    g
}

param_plot_ma <- function(
        df,
        xname,
        yname,
        mus,
        ...
    ) {

    x <- df[[xname]]
    y <- df[[yname]]

    df <- data.frame(
        M = y - x,
        mu = mus
    )
    m <- max(abs(df$M))
    ylim <- c(-m, m)

    g <- param_plot(df, "mu", "M", ...) +
        # geom_hline(yintercept = 0, lty = "twodash", col = "grey70", na.rm = TRUE) +
        ylim(ylim) +
        scale_x_log10()
    g
}

param_plot <- function(
        df,
        xname,
        yname,
        title,
        log = FALSE,
        bins = 100,
        ...
    ) {

    g <- ggplot(
        df,
        mapping = aes_string(
        x = paste_aes(xname),
        y = paste_aes(yname)
        )
    ) +
        # geom_point(
        #   alpha = 0.1, na.rm = TRUE
        # ) +
        scale_fill_viridis(direction = 1, name = "Density") +
        geom_hex(aes(fill = ..density..), bins = bins, na.rm = TRUE) +
        guides(fill = "none") +
        xlab(xname) +
        ylab(yname) +
        ggtitle(title)
    if (log) {
        g <- g +
        scale_x_log10() +
        scale_y_log10()
    }
    g
}

extract_hpd_interval_width <- function(chain, param) {
    summ <- BASiCS::Summary(chain)
    abs(summ@parameters[[param]][, 3] - summ@parameters[[param]][, 2])
}

do_hpd_plots <- function(i, maxdf, references) {
    cat(i, "/", nrow(maxdf), "\n")
    maxdf <- maxdf[i, ]
    d <- maxdf[["data"]]
    b <- maxdf[["by"]]
    rc <- references[[which(references$data == d)[[1]], "chain"]]
    cc <- readRDS(maxdf[["file"]][[1]])
    if (is.list(cc)) {
        suppressMessages(
        cc <- BASiCS:::.combine_subposteriors(
            cc,
            GeneOrder = rownames(rc),
            CellOrder = colnames(rc),
            SubsetBy = "gene"
        )
        )
    }
    cc <- BASiCS:::.offset_correct(cc, rc)
    l2 <- if (b == "advi") "ADVI" else "D & C"
    l3 <- "log2(95% HPD interval width ratio)"
    l <- list(
        plot_hpd_interval_width(
            rc,
            "Reference",
            cc,
            l2,
            "mu",
            type = "ma",
            log = FALSE,
            bins = 50
        ) +
        labs(
            title = NULL,
            x = bquote(mu[i]),
            y = l3
        ) +
        geom_hline(aes(yintercept = 0)),
        plot_hpd_interval_width(
            rc,
            "Reference",
            cc,
            l2,
            "epsilon",
            type = "ma",
            log = FALSE,
            bins = 50
        ) +
        labs(
            title = NULL,
            x = bquote(mu[i]),
            y = l3
        ) +
        geom_hline(aes(yintercept = 0))
    )
    ggsave(l[[1]],
        file = here(sprintf("figs/hpd/mu_%s_%s.pdf", b, d)),
        width = 7, height = 5
    )
    ggsave(l[[2]],
        file = here(sprintf("figs/hpd/epsilon_%s_%s.pdf", b, d)),
        width = 7, height = 5
    )
    l
}

do_ess_plot <- function(i, maxdf, references) {
    cat(i, "/", nrow(maxdf), "\n")
    maxdf <- maxdf[i, ]
    d <- maxdf[["data"]]
    b <- maxdf[["by"]]
    rc <- references[[which(references$data == d)[[1]], "chain"]]
    cc <- readRDS(maxdf[["file"]][[1]])
    if (is.list(cc)) {
        suppressMessages(
        cc <- BASiCS:::.combine_subposteriors(
            cc,
            GeneOrder = rownames(rc),
            CellOrder = colnames(rc),
            SubsetBy = "gene"
        )
        )
    }
    cc <- BASiCS:::.offset_correct(cc, rc)
    ess <- BASiCS_EffectiveSize(cc, "epsilon")

    df <- data.frame(
        x = colMedians(rc@parameters[["mu"]]),
        y = colMedians(rc@parameters[["epsilon"]] - cc@parameters[["epsilon"]]),
        color = ess
    )
    df <- df[order(df$color, decreasing = TRUE), ]
    ggplot(df, aes(x = x, y = y, color = color)) + 
        geom_point(alpha = 0.75) +
        scale_x_log10() +
        labs(x = bquote(mu[i]), y = bquote(epsilon[i]^Ref - epsilon[i]^Scalable)) +
        scale_color_viridis(
        name = bquote('Effective sample size'~epsilon[i]),
        trans = "log10"
        )
    ggsave(
        file = here(sprintf("figs/ess/%s_%s.pdf", b, d)),
        width = 5, height = 4
    )
}

paste_aes <- function(...) {
    paste0("`", ..., "`")
}


plot_hpds <- function(
        summary1,
        summary2,
        param = "mu",
        ord,
        scalename = "Normalisation    ",
        labels = c("Inferred", "Fixed")
    ) {

    df1 <- as.data.frame(summary1@parameters[[param]])
    df2 <- as.data.frame(summary2@parameters[[param]])
    # ord <- order(df1$median)
    df1 <- df1[ord, ]
    df2 <- df2[ord, ]
    df1$index <- 1:nrow(df1)
    df2$index <- 1:nrow(df2)

    if (param %in% c("mu", "delta")) {
        scale <- scale_y_log10(name = param)
    } else {
        scale <- scale_y_continuous(name = param)
    }

    g <- ggplot() +
        geom_point(
            data = df1,
            alpha = 0.6,
            size = 0.5,
            aes(x = index, y = median, colour = labels[[1]])
        ) +
        geom_point(
            data = df2,
            alpha = 0.6,
            size = 0.5,
            aes(x = index, y = median, colour = labels[[2]])
        ) +
        # geom_ribbon(
        #     data = df1,
        #     alpha = 0.4,
        #     aes(x = index, ymin = lower, ymax = upper, fill = labels[[1]])
        # ) +
        # geom_ribbon(
        #     data = df2,
        #     alpha = 0.4,
        #     aes(x = index, ymin = lower, ymax = upper, fill = labels[[2]])
        # ) +
        scale +
        xlab("Gene") +
        # guides(colour = "none") +
        theme(legend.position = "bottom") +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank()
        ) +
        scale_colour_manual(
            name = scalename,
            values = setNames(c("firebrick", "dodgerblue"), labels),
            aesthetics = c("fill", "colour")
        )
    rasterise(g, dpi = 300, layers = c("Point", "Ribbon"))
}


plot_hpd_diff <- function(
        summary_ref, summary_dc, param = "mu", ord,
        ylab = sprintf("Difference in posterior quantities\n(%s)", param)
    ) {
    df_var <- as.data.frame(summary_ref@parameters[[param]])
    df_fix <- as.data.frame(summary_dc@parameters[[param]])
    # ord <- order(df_var$median)
    df_var <- df_var[ord, ]
    df_fix <- df_fix[ord, ]
    df_diff <- df_var - df_fix
    df_diff$index <- 1:nrow(df_var)
    r <- range(c(df_diff$lower, df_diff$upper), na.rm = TRUE)
    if (max(abs(r)) > 20) {
        lims <- c(-10, 10)
    } else {
        lims <- r
    }
    g <- ggplot() +
        geom_ribbon(
            data = df_diff,
            colour = "darkgoldenrod",
            alpha = 0.5,
            aes(x = index, ymin = lower, ymax = upper)
        ) +
        geom_point(
            data = df_diff,
            alpha = 0.8,
            size = 0.5,
            aes(x = index, y = median)
        ) +
        labs(
          x = "Gene",
          y = ylab
        ) +
        coord_cartesian(ylim = lims) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank()
        )

        # +
        # guides(colour = "none") +
        # theme(legend.position = "bottom") +
        # scale_colour_manual(
        #     name = "Normalisation    ",
        #     values = c("firebrick", "dodgerblue"),
        #     aesthetics = c("fill", "colour")
        # )
    rasterise(g, dpi = 300, layers = c("Point", "Ribbon"))
}

plot_with_legend_below <- function(
        ...,
        rel_heights = c(0.85, 0.15),
        nrow = 1,
        labels = "AUTO",
        align = "h"
    ) {
    legends <- lapply(list(...), cowplot::get_legend)
    # if (!all(sapply(legends, function(x) identical(x$grobs, legends[[1]]$grobs)))) {
    #   stop("Different legends shouldn't be merged!")
    # }
    t <- theme(legend.position = "none")
    together <- cowplot::plot_grid(
        plotlist = lapply(list(...), function(x) x + t),
        labels = labels,
        nrow = nrow,
        align = align
    )
    cowplot::plot_grid(
        together, legends[[1]], nrow = 2, rel_heights = rel_heights
    )
}
