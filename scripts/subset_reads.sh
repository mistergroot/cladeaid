#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --nodelist=cbsubscb16
#SBATCH --job-name=subset
#SBATCH --output=subset.out.%j
#SBATCH --error=subset.err.%j

source ~/.bashrc
conda activate ngsLCA

INFILE1=$1
INFILE2=$2
OUTFILE1=$3
OUTFILE2=$4
NBASES=$5

cd $SLURM_SUBMIT_DIR

reformat.sh -Xmx20g \
	in1=${INFILE1} \
	in2=${INFILE2} \
	out1=${OUTFILE1} \
	out2=${OUTFILE2} \
	samplebasestarget=${NBASES} \
	sampleseed=123456
