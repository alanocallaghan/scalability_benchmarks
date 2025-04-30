conda env create -n scalability -f scalability.yml

## also need to devtools::install_github some stuff
conda activate scalability

Rscript -e 'BiocManager::install(c("BASiCS", "BiocParallel", "scater", "scran", "SingleCellExperiment", "scRNAseq"), version=3.14)'

Rscript -e 'devtools::install_github("catavallejos/BASiCS", ref="49186bda2c5e410a87c9dd8e19d20310c97e5465")'
Rscript -e 'devtools::install_github("Alanocallaghan/BASiCStan", ref="9e632610cf463c51d3856d763a89f555dc3c114c")'
Rscript -e 'devtools::install_github("jorainer/ensembldb", ref="058be8a")'
Rscript -e 'devtools::install_github("Bioconductor/BiocFileCache", ref="004cb8e")'
Rscript -e 'install.packages("patchwork")'

conda env create -n scalability-datasets -f data-env.yml
conda env create -n scalability-plotting -f plotting-env.yml
