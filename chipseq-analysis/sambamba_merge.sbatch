#!/bin/bash

## Script for running sambamba merge
## Date: 12 March 2020 

## Example usage:
## bash sambamba_merge.sh

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J sambamba_markdup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load modules
module load sambamba samtools

# Make temporary directories
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${inDir}/samtools_tmp/${SLURM_JOB_ID}

# Merge technical replicates with sambamba merge
# Might need to modify to suit your needs (i.e. four inFiles vs two)
echo $(date +"[%b %d %H:%M:%S] Running sambamba merge...")

echo "Processing sample: "${outFile/.merged.bam/}

sambamba merge \
--nthreads=${numThreads} \
${outDir}/${outFile} \
${inDir}/${inFile1} \
${inDir}/${inFile2}

# Sorting merged bam files
echo $(date +"[%b %d %H:%M:%S] Sorting merged bam...")

samtools sort \
-@ ${numThreads} \
-T ${inDir}/samtools_tmp/${SLURM_JOB_ID} \
${outDir}/${outFile} \
> ${outDir}/${outFile/.merged.bam/.merged.sorted.bam}

# Index merged bam files
echo $(date +"[%b %d %H:%M:%S] Indexing merged bam...")

samtools index \
${outDir}/${outFile/.merged.bam/.merged.sorted.bam} \
> ${outDir}/${outFile/.merged.sorted.bam/.merged.sorted.bam.bai}

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-nthreads=<int>' - number of threads to use