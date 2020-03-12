#!/bin/bash

## Script for running bwa
## Date: 31 Jan 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/trimmomatic \
## outDir=/Users/coke6162/project20xx_chipseq/bwa \
## bwaIndexDir=/scratch/Shares/chuong/genomes/hg38 \
## bwaIndex=hg38.main.fa \
## sbatch --array 0-42 bwa_PE.sbatch

## bwaIndexDir options: hg38, mm10, bosTau8, canFam3
## bwaIndex options: hg38.main.fa, mm10.fa, bosTau8.fa, canFam3.fa 

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=8:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load bwa and samtools
module load bwa/0.7.15 samtools/1.8

# Define query files
# Removes R1/R2.fastq.gz to generate a unique identifier for each pair of files
# Modify to fit input file name if necessary
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename | sed --expression='s/_R1_bbduk.fastq.gz//g' | sed --expression='s/_R2_bbduk.fastq.gz//g' | uniq))

# Make temporary directories
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${outDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Run bwa
echo "Aligning reads to the genome & outputting sam files"
echo "bwa version: "$(bwa)
echo "Samtools version: "$(samtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bwa mem \
-t ${numThreads} \
${bwaIndexDir}/${bwaIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_bbduk.fastq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_bbduk.fastq.gz \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

# Generate flagstat files from unfiltered sam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from unfiltered sam...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_flagstat.txt

# Convert sam files to bam files, filtering out low quality (<10) and unmapped (through '-F 4' option)
echo $(date +"[%b %d %H:%M:%S] Converting sam to unfiltered bam & sorting")

samtools view \
-@ ${numThreads} \
-Sb \
-q 10 \
-F 4 \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
| samtools sort \
-@ ${numThreads} \
-T ${outDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} \
- \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Index sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing sorted bam...")

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing sam...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-t <int' - specifies number of threads to use
# '-Sb' - convert sam (bwa output) to bam
# '-q <int>' - skip alignments with MAPQ (i.e. quality) score smaller than <int>
# '-@ <int>' - specifies number of threads to use