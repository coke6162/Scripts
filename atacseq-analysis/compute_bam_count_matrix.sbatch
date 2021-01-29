#!/bin/bash

## Script for running bowtie2 index
## Date: 29 January 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## bedFile=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/DESeq2/exampleBed.bed \
## refFile=/Users/coke6162/repos/project_files/yyyymmdd_ATAC-seq_projectName/samples.txt \
## bamDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/bowtie2/dedup/shifted \
## outDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/DESeq2 \
## name=narrowPeak \
## sbatch compute_bam_count_matrix.sbatch

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

# create a file containing full paths to all input bams
pwd; hostname; date

echo "Creating a file containing full paths to all input bams..."

for i in $(cat ${refFile})
do 
	ls ${bamDir}/${i}*.bam
done > ${outDir}/input_bams.txt

# Compute a matrix of bam alignment counts for each region of interest
pwd; hostname; date

echo "Starting bedtools..."

echo $(date +"[%b %d %H:%M:%S] Computing matrix of bam alignment counts...")

bedtools multicov \
-bams $(cat ${outDir}/input_bams.txt) \
-bed ${bedFile} \
> ${outDir}/bam_counts_${name}_d100.tab

echo $(date +"[%b %d %H:%M:%S] Done")