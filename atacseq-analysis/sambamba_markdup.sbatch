#!/bin/bash

## Script for running sambamba markdup
## Date: 1 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## inDir=/Users/coke6162/project20xx_atacseq/bowtie2 \
## outDir=/Users/coke6162/project20xx_atacseq/bowtie2/dedup \
## refFile=/Users/CL_Shared/db/genomes/canFam3/fa/canFam3.fa \
## sbatch --array=0-3 sambamba_markdup.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=16G

# Job name and output
#SBATCH -J sambamba_markdup
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
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename))

# Filter out reads mapping to mitochondrial genes
pwd; hostname; date

echo "Removing reads mapping to mitochondrial genes"
echo "Samtools version: "$(samtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

samtools view \
-@ ${numThreads} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| egrep -v chrM \
| samtools view -bT ${refFile} -o ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam

echo $(date +"[%b %d %H:%M:%S]")

# Sort bam files and remove duplicates with sambamba
echo "Sorting bam files & removing duplicates with sambamba"
echo "Sambamba version: see stdout"`sambamba --version`

sambamba sort \
--natural-sort \
--nthreads=${numThreads} \
--memory-limit=4G \
--tmpdir=${inDir}/tmp \
--out=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam

sambamba markdup \
--remove-duplicates \
--nthreads=${numThreads} \
--overflow-list-size=600000 \
--tmpdir=${inDir}/tmp \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.dedup.bam

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam

echo $(date +"[%b %d %H:%M:%S]")

# Run samtools flagstat 
echo "Running flagstat"

samtools flagstat \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.dedup.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.dedup.flagstat.txt

echo $(date +"[%b %d %H:%M:%S]")

# Resort dedup output
echo "Sort dedup bams"

samtools sort \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.dedup.bam \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.dedup.resorted.bam

# Index resorted bams
echo "Index resorted bams"

samtools index \
-@ {numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.dedup.resorted.bam

echo $(date +"[%b %d %H:%M:%S] Done")