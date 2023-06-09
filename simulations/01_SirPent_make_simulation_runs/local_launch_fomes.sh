#! /bin/bash

ROOT=/Users/nbrazeau/Documents/Github/PFSwDemonDanProphetBrown/ # root directory for project (non-scratch)
WD=/Users/nbrazeau/Documents/Github/PFSwDemonDanProphetBrown/data/temp # scratch directory for networks
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/simulations/01_SirPent_make_simulation_runs/run_fomes.snake \
	--configfile $ROOT/simulations/01_SirPent_make_simulation_runs/config_fomes_batch_local.yaml \
	--printshellcmds \
	--directory $WD \
	--cores 1 \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--dryrun -p
