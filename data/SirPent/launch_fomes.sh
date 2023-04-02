#! /bin/bash

ROOT=*** # root directory for project (non-scratch)
WD=*** # scratch directory for networks
NODES=1028 # max number of cluster nodes
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/run_fomes.snake \
	--configfile config_fomes_batch.yaml \
	--printshellcmds \
	--directory $WD \
	--cluster $ROOT/launch.py \
	-j $NODES \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT
