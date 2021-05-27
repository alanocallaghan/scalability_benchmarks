do_fit_plot <- function(i, maxdf, references) {
  c <- readRDS(maxdf[[i, "file"]])
  b <- maxdf[[i, "by"]]
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        SubsetBy = "gene"
      )
    )
  }
  d <- maxdf[[i, "data"]]
  rc <- references[[which(references$data == d), "chain"]]
  l2 <- if (b == "advi") "ADVI" else "Divide and conquer"
  c <- BASiCS:::.offset_correct(rc, c)
  g1 <- BASiCS_ShowFit(rc) 
  g2 <- BASiCS_ShowFit(c) 
  ggsave(g1, 
    file = here(paste0("figs/fit/", d, "_ref", ".pdf")),
    width = 4,
    height = 3
  )
  ggsave(g2,
    file = here(paste0("figs/fit/", d, "_", b, ".pdf")), 
    width = 4,
    height = 3
  )
  cowplot::plot_grid(
    g1 + ggtitle("Reference"),
    g2 + ggtitle(l2)
  )
}

do_de_plot <- function(i, maxdf, references) {
  d <- maxdf[[i, "data"]]
  b <- maxdf[[i, "by"]]
  rc <- references[[which(references$data == d), "chain"]]

  c <- readRDS(maxdf[[i, "file"]])
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        GeneOrder = rownames(rc),
        CellOrder = colnames(rc),
        SubsetBy = "gene"
      )
    )
  }
  l2 <- if (b == "advi") "ADVI" else "D & C"
  de <- BASiCS_TestDE(rc, c,
    GroupLabel1 = "Ref",
    GroupLabel2 = l2,
    Plot = FALSE,
    PlotOffset = FALSE,
    EFDR_M = NULL,
    EFDR_D = NULL,
    EFDR_R = NULL
  )
  g <- BASiCS_PlotDE(de@Results[[1]],
    Plots = c("MAPlot", "VolcanoPlot"),
    Mu = de@Results$Mean@Table$MeanOverall
  )
  ggsave(g, file = here(paste0("figs/de/mu_", d, "_", b, ".pdf")), width = 9, height = 5)
  g <- BASiCS_PlotDE(
    de@Results[[3]],
    Plots = c("MAPlot", "VolcanoPlot"),
    Mu = de@Results$Mean@Table$MeanOverall
  )
  ggsave(g, file = here(paste0("figs/de/epsilon_", d, "_", b, ".pdf")), width = 9, height = 5)
  g
}

plot_hpd_interval <- function(
    chain1,
    xname,
    chain2,
    yname,
    param,
    type = c("ma", "std"),
    log = FALSE,
    ...
  ) {

  x <- extract_hpd_interval(chain1, param)
  y <- extract_hpd_interval(chain2, param)

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
      chain1@parameters[["mu"]][, 1],
      chain2@parameters[["mu"]][, 1]
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
      "of HPD interval\nof", param, "parameter in",
      yname, "vs", xname, "MCMC"
    ),
    log = log && type =="std",
    mus = mus,
    ...
  )

  # ggMarginal(g)
  g
}

extract_hpd_interval <- function(chain, param) {
  summ <- BASiCS::Summary(chain)
  abs(summ@parameters[[param]][, 3] - summ@parameters[[param]][, 2])
}

do_hpd_plots <- function(j, maxdf, references) {
  d <- maxdf[[j, "data"]]
  b <- maxdf[[j, "by"]]
  rc <- references[[which(references$data == d), "chain"]]
  c <- readRDS(maxdf[[j, "file"]])
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        GeneOrder = rownames(rc),
        CellOrder = colnames(rc),
        SubsetBy = "gene"
      )
    )
  }
  c <- BASiCS:::.offset_correct(rc, c)
  l2 <- if (b == "advi") "ADVI" else "D & C"
  l3 <- paste0("log2(", l2, " HPD interval width / Reference HPD interval width)")
  l <- list(
    plot_hpd_interval(
      rc,
      "Reference",
      c,
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
    plot_hpd_interval(
      rc,
      "Reference",
      c,
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
  ggsave(l[[1]], file = here(paste0("figs/hpd/mu_", b, "_", d, ".pdf")), width = 7, height = 5)
  ggsave(l[[2]], file = here(paste0("figs/hpd/epsilon_", b, "_", d, ".pdf")), width = 7, height = 5)
  l
}

do_ess_plot <- function(j, maxdf, references) {
  d <- maxdf[[j, "data"]]
  b <- maxdf[[j, "by"]]
  rc <- references[[which(references$data == d), "chain"]]
  c <- readRDS(maxdf[[j, "file"]])
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        GeneOrder = rownames(rc),
        CellOrder = colnames(rc),
        SubsetBy = "gene"
      )
  }
  c <- BASiCS:::.offset_correct(rc, c)

  ess <- effectiveSize(as.mcmc(c@parameters[["epsilon"]]))

  df <- data.frame(
    x = colMedians(rc@parameters[["mu"]]),
    y = colMedians(rc@parameters[["epsilon"]] - c@parameters[["epsilon"]]),
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
  ggsave(file = here(paste0("figs/ess/", b, "_", d, ".pdf")), width = 5, height = 4)
}





dir.create("figs/ess", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/hpd", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/fit", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/de", showWarnings = FALSE, recursive = TRUE)



## Mean-variance curves and DE plots for worst for each data
maxdfm <- df %>% filter(!is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffExp) %>%
  distinct(data, .keep_all = TRUE)


maxdfe <- df %>% filter(!is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffResDisp) %>%
  distinct(data, .keep_all = TRUE)


ess_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_ess_plot,
  maxdf = maxdfe,
  references = references
)

fit_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_fit_plot,
  maxdf = maxdfe,
  references = references
)

## MA plot of worst case for worst for each data
de_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_de_plot,
  maxdf = maxdfe,
  references = references
)

hpd_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_hpd_plots,
  maxdf = maxdfe,
  references = references
)

## Posterior variance plots
maxdfadvim <- df %>% filter(is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffExp) %>%
  distinct(data, .keep_all = TRUE)

maxdfadvie <- df %>% filter(is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffResDisp) %>%
  distinct(data, .keep_all = TRUE)


fit_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_fit_plot,
  maxdf = maxdfadvie, 
  references = references
)

## MA plot of worst case for worst for each data
de_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_de_plot,
  maxdf = maxdfadvie,
  references = references
)

hpd_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_hpd_plots,
  maxdf = maxdfadvie,
  references = references
)