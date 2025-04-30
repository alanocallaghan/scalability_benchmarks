# grep -rh "library(" | sed -e 's/^[ \t]*//' | sort | uniq
# library("argparse")
# library("coda")
# library("dplyr")
# library("ggbeeswarm")
# library("ggplot2")
# library("here")
# library("viridis")
# library("BASiCStan")
# library("BASiCS")
# library("BiocParallel")
# library("scater")
# library("scran")
# library("SingleCellExperiment")


conda create -n scalability \
    gcc \
    r-base=4.1.1 \
    r-argparse \
    r-curl \
    r-httr \
    r-coda \
    r-dplyr \
    r-ggbeeswarm \
    r-ggplot2 \
    r-here \
    r-viridis \
    snakemake \
    r-biocmanager \
    r-devtools \
    r-rstan \
    r-ggpointdensity \
    r-rcpparmadillo \
    r-ggrastr \
    r-patchwork \
    bioconductor-scrnaseq
    # bioconductor-basics \
    # bioconductor-scater \
    # bioconductor-scran \
    # bioconductor-singlecellexperiment \


## also need to devtools::install_github some stuff

# mamba activate scalability
Rscript -e 'BiocManager::install(c("BASiCS", "BiocParallel", "scater", "scran", "SingleCellExperiment", "scRNAseq"), version=3.14)'

Rscript -e 'devtools::install_github("catavallejos/BASiCS", ref="49186bda2c5e410a87c9dd8e19d20310c97e5465")'
Rscript -e 'devtools::install_github("Alanocallaghan/BASiCStan", ref="9e632610cf463c51d3856d763a89f555dc3c114c")'
Rscript -e 'devtools::install_github("jorainer/ensembldb", ref="058be8a")'
Rscript -e 'devtools::install_github("Bioconductor/BiocFileCache", ref="004cb8e")'
Rscript -e 'install.packages("patchwork")'
