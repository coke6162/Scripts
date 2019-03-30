#!/bin/bash

## Script for running hisat2
## Date: 30 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_rnaseq/trimmomatic \
## outDir=/Users/coke6162/project20xx_rnaseq/hisat2 \
## hisatIndexDir=/scratch/Shares/chuong/genomes/hg38 \
## hisatIndex=hg38.main \
## sbatch --array 0-42 hisat2_SE.sbatch

## hisatIndexDir options: hg38, mm10
## hisatIndex: hg38.main, mm10.main

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J hisat2_SE
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@colorado.edu

# Set constant variables
numThreads=8

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls ${inDir}/*_trimmed.fastq.gz | xargs -n 1 basename))

# Load hisat2 and samtools
module load hisat2 samtools

# Run hisat2
pwd; hostname; date

echo "HISAT2 version: "$(hisat2 --version)
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Starting HISAT2 alignment..."

hisat2 -p ${numThreads} --no-softclip -x ${hisatIndexDir}/${hisatIndex} \
-U ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_trimmed.fastq.gz}.sorted.uniq.bam

# Create .sorted.uniq.bai index files for upload to UCSC browser
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_trimmed.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-p <int>' - number of threads
# '--no-softclip' - disallows soft clipping
# '-x <file.basename>' - hisat2 index basename; include everything up to .main.1.ht2
# '-U <file>' - input single-end data
# '-q <int>' - filter out mapped reads by quality threshold
# '-Sb' - convert sam to bam
# '-@' - number of threads