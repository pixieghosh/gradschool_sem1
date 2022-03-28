#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --time=5:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=samtools


module purge
module load samtools/intel/1.11 

Name=$1

#samtools view -b $Name.sam > $Name.bam
#samtools sort $Name.bam -o $Name.sorted.bam
samtools index $Name.sorted.bam