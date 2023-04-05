#! /bin/bash

ROOT=/Users/nbrazeau/Documents/Github/PFSwDemonDanProphetBrown/ # root directory for project (non-scratch)
WD=/Users/nbrazeau/Documents/Github/PFSwDemonDanProphetBrown/data/temp # scratch directory for networks
NODES=1028 # max number of cluster nodes
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/data/SirPent/run_fomes.snake \
	--configfile $ROOT/data/SirPent/config_fomes_batch.yaml \
	--printshellcmds \
	--directory $WD \
	--cluster $ROOT/data/SirPent/launch.py \
	-j $NODES \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
#	--dryrun -p
