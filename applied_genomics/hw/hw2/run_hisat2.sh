#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=hisat


module purge
module load hisat2/2.2.1

index=$1
LREAD=$2
RREAD=$3

hisat2-build -f $index.fa $index 

hisat2 --threads=6 -x $index -1 $LREAD -2 $RREAD > $LREAD.sam