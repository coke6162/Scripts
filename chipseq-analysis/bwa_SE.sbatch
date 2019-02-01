#!/bin/bash

## Script for running bwa
## Date: 31 Jan 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/trimmomatic \
## outDir=/Users/coke6162/project20xx_chipseq/bwa \
## bwaIndexDir=/scratch/Shares/chuong/genomes/hg38 \
## bwaIndex=hg38.main.fa \
## sbatch --array 0-42 bwa_SE.sbatch

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
queries=($(ls ${inDir}/*_trimmed.fastq.gz | xargs -n 1 basename))

# Run bwa
pwd; hostname; date

echo "Starting bwa alignment..."
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bwa mem -t 8 ${bwaIndexDir}/${bwaIndex} ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools view -Sb -q 10 - \
| samtools sort -@ 8 - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_trimmed.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-t <int' - number of threads to use
# '-Sb' - convert sam (bwa output) to bam
# '-q <int>' - skip alignments with MAPQ (i.e. quality) score smaller than <int>
# '-@ <int>' - number of threads to use