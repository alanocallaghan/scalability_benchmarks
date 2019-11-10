#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=5G
#$ -pe sharedmem 16
#$ -cwd
#$ -l h_rt=24:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.1

settings=($(sed -n "$SGE_TASK_ID p" data/downsampling_ref_grid.txt))

Rscript scripts/chain-scripts/downsampling_reference.R ${settings[@]} outputs/downsampling/reference/$SGE_TASK_ID