#! /bin/bash

ROOT=/proj/emchlab/users/NickB/Projects/PFSwDemonDanProphetBrown # root directory for project (non-scratch)
WD=/work/users/n/f/nfb/Projects/PFSwDemonDanProphetBrown/ # scratch directory for networks
NODES=1028 # max number of cluster nodes
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/simulations/01_SirPent_make_simulation_runs/run_fomes.snake \
	--configfile $ROOT/simulations/01_SirPent_make_simulation_runs/config_fomes_batch_LL.yaml \
	--printshellcmds \
	--directory $WD \
	--cluster $ROOT/simulations/01_SirPent_make_simulation_runs/launch_heavy.py \
	-j $NODES \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
#	--dryrun -p
