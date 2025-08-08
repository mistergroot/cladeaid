#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=regular
#SBATCH --nodelist=cbsubscb16
#SBATCH --job-name=trimga
#SBATCH --output=trimga.out.%j
#SBATCH --error=trimga.err.%j

source ~/.bashrc
conda activate aquaculture

INFILE1=$1
INFILE2=$2
OUTPATH=$3
OUTPREFIX=$4

cd $SLURM_SUBMIT_DIR
mkdir -p ${OUTPATH}

date
trim_galore --j ${SLURM_NTASKS} --gzip \
	--output_dir ${OUTPATH} \
	--basename ${OUTPREFIX} \
	--paired ${INFILE1} ${INFILE2}
date
