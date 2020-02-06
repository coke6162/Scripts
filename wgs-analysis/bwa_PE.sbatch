#!/bin/bash

## Script for running bwa
## Date: 6 February 2020 

## Example usage:
## cd /Users/coke6162/repos/project_files/20200102_DNA-seq_bosTau9_USMARC
## bash runBWA.sh

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=48:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load bwa
module load samtools

# Set constant variables
numThreads=8

# Define path variables
inDir=${projectDir}/fastq
outDir=${projectDir}/bwa

# Run bwa
pwd; hostname; date

echo "bwa version: "$(bwa)
echo "Samtools version: "$(samtools --version)

echo "Processing sample: "${inDir}/${fileName}

echo "Aligning reads to the genome & outputting sam files"

bwa mem \
-M \
-t ${numThreads} \
-R "@RG\tID:${fileName}\tLB:IlluminaExome\tPL:ILLUMINA\tSM:${sample}" \
${bwaIndex} \
${inDir}/${fileName}_1.fastq.gz \
${inDir}/${fileName}_2.fastq.gz \
> ${outDir}/${fileName}.sam

# Generate flagstat files from unfiltered sam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from sam...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/${fileName}.sam \
> ${outDir}/flagstat/${fileName}_flagstat.txt

# Convert sam files to bam files & sort
echo $(date +"[%b %d %H:%M:%S] Converting sam to bam & sorting")

samtools view \
-@ ${numThreads} \
-Sb \
${outDir}/${fileName}.sam \
| samtools sort \
-@ ${numThreads} \
- \
> ${outDir}/${fileName}.sorted.bam

# Index sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing sorted bams...")

samtools index \
${outDir}/${fileName}.sorted.bam \
> ${outDir}/${fileName}.sorted.bam.bai

# Remove intermediate files
rm ${outDir}/${fileName}.sam

echo $(date +"[%b %d %H:%M:%S] Done")