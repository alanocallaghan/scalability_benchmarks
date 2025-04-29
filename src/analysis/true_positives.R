library("BASiCS")
library("ggplot2")
library("here")

source(here("src/analysis/preamble.R"))
source(here("src/analysis/functions.R"))

files <- list.files(
    "outputs/true-positives/",
    pattern = ".*.rds",
    full.names = TRUE
)
files_advi <- files[grepl("advi", files)]
files <- files[!grepl("advi", files)]
# todo: advi bit

pos_metadata <- data.frame(
    data = "ibarra-soria",
    nsubsets = gsub(".*nsubsets-(\\d+).*", "\\1", files),
    seed = gsub(".*seed-(\\d+).*", "\\1", files),
    file = files
)

pos_metadata_advi <- data.frame(
    data = "ibarra-soria",
    nsubsets = 0,
    seed = gsub(".*-(\\d+).rds", "\\1", files_advi),
    file = files_advi
)
pos_metadata <- rbind(pos_metadata, pos_metadata_advi)

pos_metadata$test <- parallel::mclapply(seq_along(pos_metadata$file),
    function(i) {
        cat(i, "/", nrow(pos_metadata), "\n")
        readRDS(pos_metadata$file[[i]])$test
    }, mc.cores = 12
)

params <- c("Mean", "Disp", "ResDisp")
pg <- paste0(params, "Gene")
pos_metadata[, pg] <- NA

pos_metadata[, pg] <- lapply(
    pos_metadata[, pg],
    function(x) rep(list(), length(x))
)


for (i in seq_len(nrow(pos_metadata))) {
    cat(i, "/", nrow(pos_metadata), "\n")
    df <- lapply(
        params,
        function(x) {
            out <- as.data.frame(pos_metadata$test[[i]], Parameter = x)$GeneName
            if (length(out)) {
                out
            } else {
                NA
            }
        }
    )
    for (j in seq_along(pg)) {
        pos_metadata[[i, pg[[j]]]] <- list(df[[j]])
    }
}
pos_metadata$test <- NULL

ref <- pos_metadata[which(pos_metadata$nsubsets == 1), ]
pos_metadata <- pos_metadata[pos_metadata$nsubsets != 1, ]
jp <- paste0(params, "Jaccard")
pos_metadata[, jp] <- NA
for (i in seq_len(nrow(pos_metadata))) {
    cat(i, "/", nrow(pos_metadata), "\n")
    for (j in seq_along(jp)) {
        pos_metadata[i, jp[[j]]] <- jaccard(
            ref[[1, pg[[j]]]][[1]],
            pos_metadata[[i, pg[[j]]]][[1]]
        )
    }
}



jd <- pos_metadata[, c("nsubsets", jp)]
mdf_tp <- reshape2::melt(jd, id.var = "nsubsets")
mdf_tp$variable <- gsub("Jaccard", "", mdf_tp$variable)
mdf_tp$nsubsets <- factor(
    mdf_tp$nsubsets,
    levels = sort(unique(as.numeric(mdf_tp$nsubsets)))
)

mdf_tp$variable <- plyr::revalue(
    mdf_tp$variable,
    replace = c(
        "Mean" = "mu",
        "Disp" = "delta",
        "ResDisp" = "epsilon"
    )
)
mdf_tp$variable <- factor(
    mdf_tp$variable,
    levels = c("mu", "delta", "epsilon")
)
mdf_tp <- mdf_tp[mdf_tp$variable %in% c("mu", "epsilon"), ]
mdf_tp$variable <- droplevels(mdf_tp$variable)
levels(mdf_tp$variable) <- c("Mean", "Residual over-dispersion")


## Mean-variance curves and DE plots for worst for each data
# maxdfm <- mdf_tp %>%
#   group_by(data) %>%
#   top_n(n = 1, wt = -value) %>%
#   distinct(data, .keep_all = TRUE)


mdf_tp_advi <- mdf_tp %>%
    filter(nsubsets == 0) %>%
    group_by(variable) %>%
    summarise(value = mean(value))
mdf_tp <- mdf_tp[mdf_tp$nsubsets != 0, ]

g <- ggplot() +
    geom_quasirandom(
        data = mdf_tp,
        aes(nsubsets, y = value, colour = variable),
        dodge.width = 0.3, size = 0.7, width = 0.1, groupOnX = TRUE
    ) +
    geom_hline(
        data = mdf_tp_advi,
        linetype = "dashed",
        aes(colour = variable, yintercept = value),
        show_guide=FALSE,
    ) +
    labs(x = "Number of subsets", y = "Jaccard Index") +
    ylim(0, 1) +
    scale_color_brewer(palette = "Dark2", name = "Parameter") +
    theme_bw() +
    theme(legend.position = "bottom") +
    theme(panel.grid = element_blank())

ggsave(file = "figs/true_positives.pdf", width = 5, height = 3.5)
save.image("rdata/true_positives.RData")
