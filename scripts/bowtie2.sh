#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=regular
#SBATCH --nodelist=cbsubscb16
#SBATCH --job-name=bowtiemap
#SBATCH --output=bowtie.out.%j
#SBATCH --error=bowtie.err.%j

source ~/.bashrc
conda activate ngsLCA

INFILE1=$1
INFILE2=$2
REFERENCE=$3
OUTPREFIX=$4

cd $SLURM_SUBMIT_DIR
mkdir -p ${OUTPATH}

bowtie2 --threads ${SLURM_NTASKS} -k 500 --end-to-end --no-unal \
    -x ${REFERENCE} \
    -1 ${INFILE1} -2 ${INFILE2} \
    | awk '$3!="*"' \
    | samtools view -bS - > ${OUTPREFIX}.bam
