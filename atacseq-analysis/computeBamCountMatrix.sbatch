#!/bin/bash

## Script for running bowtie2 index
## Date: 17 June 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## bedFile=/Users/coke6162/project20xx_atacseq/DESeq2/exampleBed.bed \
## bamFile=/Users/coke6162/repos/project_files/project20xx_atacseq/exampleBam.bam \
## name=narrowPeak \
## outDir=/Users/coke6162/project20xx_atacseq/DESeq2 \
## sbatch computeBamCountMatrix.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J computeBamCountMatrix
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load bedtools/2.25.0

# define key variables
bams=$(cat ${bamFile})

# Compute a matrix of bam alignment counts for each region of interest
pwd; hostname; date

echo "Starting bedtools..."

echo $(date +"[%b %d %H:%M:%S] Computing matrix of bam alignment counts...")

bedtools multicov \
-bams ${bams} \
-bed ${bedFile} \
> ${outDir}/bamCounts_${name}_d100.tab

echo $(date +"[%b %d %H:%M:%S] Done")