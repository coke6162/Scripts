#!/bin/bash

## Script for running bedtools slop
## Date: 1 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## inDir=/Users/coke6162/project20xx_atacseq/macs2 \
## outDir=/Users/coke6162/project20xx_atacseq/macs2/OCR \
## chromSizes=/Users/CL_Shared/db/genomes/bosTau8/fa/bosTau8.chrom.sizes \
## sbatch --array=0-3 bowtie2.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=4G

# Job name and output
#SBATCH -J bedtools_slop
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Module load bowtie2 and samtools
module load bedtools

# Define query files
queries=($(ls ${inDir}/*_summits.bed | xargs -n 1 basename))

# Extend summits 125bp each way
pwd; hostname; date

echo "Slopping summits to get 250bp OCRs"
echo "bedtools version: "$(bedtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bedtools slop \
-i ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-g ${chromSizes} \
-b 125 > ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_summits.bed}_250bp_OCR.bed

echo $(date +"[%b %d %H:%M:%S] Done")
