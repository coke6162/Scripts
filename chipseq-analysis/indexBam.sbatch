#!/bin/bash

## Script for indexing bams
## Date: 11 July 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/chipseq-analysis
## inDir=/Users/coke6162/project20xx_chipseq/bwa \
## outDir=/Users/coke6162/project20xx_chipseq/bwa \
## sbatch --array=0-3 indexBam.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J indexBam
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load modules
module load samtools

# Define query files
queries=($(ls ${inDir}/*.sorted.uniq.bam | xargs -n 1 basename))

# Index sorted bam files
pwd; hostname; date

echo "Indexing sorted bams"
echo "Samtools version: "$(samtools --version)

samtools index \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bai

echo $(date +"[%b %d %H:%M:%S] Done")