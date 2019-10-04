#!/bin/bash

## Script for running samtools merge
## Date: 15 July 2019 

## Example usage:
## inDir1=/Users/coke6162/project20xx_rnaseq/hisat2 \
## inDir2=/Users/coke6162/project20xx_rnaseq/hisat2 \
## outDir=/Users/coke6162/project20xx_rnaseq/hisat2/merged \
## name=test \
## bam1=test1.sorted.uniq.bam \
## bam2=test2.sorted.uniq.bam \
## sbatch mergeBams.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J mergeBams
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load modules
module load samtools

# Merge bams
pwd; hostname; date

echo "samtools version: "$(samtools --version)
echo "Merging bams: "$(echo ${bam1} ${bam2})

echo $(date +"[%b %d %H:%M:%S] Merging bams...")

samtools merge \
-@ ${numThreads} \
${outDir}/${name}_merged.sorted.uniq.bam \
${inDir1}/${bam1} \
${inDir2}/${bam2}

# Index merged bams
echo $(date +"[%b %d %H:%M:%S] Indexing merged bams")

samtools index ${outDir}/${name}_merged.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")