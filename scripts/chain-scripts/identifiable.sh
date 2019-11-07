#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -cwd
#$ -l h_rt=24:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.1

settings=($(sed -n "$SGE_TASK_ID p" data/datasets_batch.txt))

Rscript scripts/chain-scripts/identifiable.R ${settings[@]} /exports/eddie/scratch/s1372510/identifiable/$SGE_TASK_ID
Rscript scripts/chain-scripts/non_identifiable.R ${settings[@]} /exports/eddie/scratch/s1372510/identifiable/$SGE_TASK_ID
