#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=regular
#SBATCH --nodelist=cbsubscb16
#SBATCH --job-name=fastqdump
#SBATCH --output=fastqdump.out.%j
#SBATCH --error=fastqdump.err.%j

source ~/.bashrc
conda activate cladeaid

ACCESSION=$1
SPOTS=$2

cd $SLURM_SUBMIT_DIR
fastq-dump.3.2.1 \
    ${ACCESSION} \
    -X ${SPOTS} \
    --split-3
