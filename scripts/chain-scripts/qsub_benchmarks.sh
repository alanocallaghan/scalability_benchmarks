# qsub -t 1-$(cat data/divide_and_conquer_grid.txt | wc -l) -tc 10 scripts/chain-scripts/divide_and_conquer.sh
# qsub -t 1-$(cat data/advi_grid.txt | wc -l) -tc 5 scripts/chain-scripts/advi.sh
# qsub -t 1-$(cat data/downsampling_ref_grid.txt | wc -l) -tc 5 scripts/chain-scripts/downsampling_reference.sh
# qsub -t 1-$(cat data/downsampling_grid.txt | wc -l) -tc 5 scripts/chain-scripts/downsampling_divide.sh
qsub -t 1-$(cat data/datasets_batch.txt | wc -l) -tc 5 scripts/chain-scripts/identifiable.sh
qsub -t 1-$(cat data/datasets_batch.txt | wc -l) -tc 5 scripts/chain-scripts/batchinfo.sh
