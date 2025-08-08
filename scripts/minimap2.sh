#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=regular
#SBATCH --nodelist=cbsubscb16
#SBATCH --job-name=minimap
#SBATCH --output=minimap.out.%j
#SBATCH --error=minimap.err.%j

source ~/.bashrc
conda activate ngsLCA

INFILE1=$1
INFILE2=$2
REFERENCE=$3
OUTPREFIX=$4

cd $SLURM_SUBMIT_DIR
mkdir -p ${OUTPATH}

minimap2 -ax sr -t ${SLURM_NTASKS} -N 500 --sam-hit-only \
    ${REFERENCE} \
    ${INFILE1} ${INFILE2} \
    | awk '$3!="*"' \
    | samtools view -bS - > ${OUTPREFIX}.bam
