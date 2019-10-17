#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=20G
#$ -cwd
#$ -l h_rt=12:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.0

settings=($(sed -n "$SGE_TASK_ID p" data/advi_grid.txt))

Rscript scripts/advi.R ${settings[@]} /exports/eddie/scratch/s1372510/advi/$SGE_TASK_ID
