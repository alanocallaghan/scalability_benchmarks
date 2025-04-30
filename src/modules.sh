#!/bin/sh

if [ -f "/etc/profile.d/modules.sh" ]; then

    . /etc/profile.d/modules.sh

    # module load phys/compilers/gcc/10.2.0

    # export C_INCLUDE_PATH=/exports/igmm/software/pkg/el7/apps/hdf5/1.8.13/include

    # module load roslin/gcc/7.3.0

fi

set -eu
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
fi
set +eu

#conda activate $(find .snakemake/conda/ -mindepth 1 -maxdepth 1 -type d)
conda activate scalability
