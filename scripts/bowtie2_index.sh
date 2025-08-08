#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --nodelist=cbsubscb16
#SBATCH --job-name=bowtie2index
#SBATCH --output=bowtie2index.out.%j
#SBATCH --error=bowtie2index.err.%j

source ~/.bashrc
conda activate aquaculture

REFERENCE=$1

cd $SLURM_SUBMIT_DIR
bowtie2-build \
	--threads ${SLURM_NTASKS} \
	${REFERENCE} \
	${REFERENCE}
