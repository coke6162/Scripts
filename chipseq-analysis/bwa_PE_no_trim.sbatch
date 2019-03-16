#!/bin/bash

## Script for running with input that was not trimmed (i.e. didn't need to be trimmed)
## Date: 14 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/fastq \
## outDir=/Users/coke6162/project20xx_chipseq/bwa \
## bwaIndexDir=hg38 \
## bwaIndex=hg38.main.fa \
## sbatch --array 0-42 bwa_PE_no_trim.sbatch

## bwaIndexDir options: hg38, mm10, bosTau8, canFam3
## bwaIndex options: hg38.main.fa, mm10.fa, bosTau8.fa, canFam3.fa 

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=8:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load bwa and samtools
module load bwa/0.7.15 samtools/1.8

# Define query files
# Removes R1/R2.fastq.gz to generate a unique identifier for each pair of files
# Modify to fit input file name if necessary
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename | sed --expression='s/_R1.fastq.gz//g' | sed --expression='s/_R2.fastq.gz//g' | uniq))

# Run bwa
pwd; hostname; date

echo "bwa version: 0.7.15-r1140"
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Starting bwa alignment...")
echo $(date +"[%b %d %H:%M:%S] Converting bams to sams, filtering out low quality reads, and sorting...")

bwa mem -t 8 ${bwaIndexDir}/${bwaIndex} ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1.fastq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2.fastq.gz \
| samtools view -Sb -q 10 - \
| samtools sort -@ 8 - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-t <int' - specifies number of threads to use
# '-Sb' - convert sam (bwa output) to bam
# '-q <int>' - skip alignments with MAPQ (i.e. quality) score smaller than <int>
# '-@ <int>' - specifies number of threads to use