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

if [ -f "/etc/profile.d/modules.sh" ]; then
    module load roslin/gcc/7.3.0
fi

mamba create -y -n scalability2 \
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

conda activate scalability2
Rscript -e 'BiocManager::install(c("BASiCS", "BiocParallel", "scater", "scran", "SingleCellExperiment", "scRNAseq"), version=3.14)'

Rscript -e 'devtools::install_github("catavallejos/BASiCS", ref="11a03083be88d9900a168cd4ef7f41367f6fa4ed")'
Rscript -e 'devtools::install_github("Alanocallaghan/BASiCStan", ref="9e632610cf463c51d3856d763a89f555dc3c114c")'
Rscript -e 'devtools::install_github("jorainer/ensembldb")'
