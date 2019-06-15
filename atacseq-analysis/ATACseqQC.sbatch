#!/bin/bash

## Script for running ATACseqQC
## Date: 3 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## projectDir=/Users/coke6162/repos/project_files/project20xx_atacseq \
## scriptDir=/Users/coke6162/repos/scripts/atacseq-analysis \
## inDir=/Users/coke6162/project20xx_atacseq/bowtie2 \
## outDir=/Users/coke6162/project20xx_atacseq/bowtie2/ATACseqQC \
## sbatch --array=0-3 ATACseqQC.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J ATACseqQC
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define constant variables
libraryDir=/Users/coke6162/R/x86_64-redhat-linux-gnu-library/3.5/ATACseqQC/extdata

# Define query files
queries=($(ls ${inDir}/*.sam | xargs -n 1 basename | sed 's/.sam//g'))

# Load modules
module load samtools
module load R/3.5.1

# Specify R library paths
export R_LIBS='/Users/coke6162/R/x86_64-redhat-linux-gnu-library/3.5/'

# Copy bam and bai files to the ATACseqQC directory
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

echo "Copying bam and bai files to ATACseqQC directory"

cp ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam ${libraryDir}
cp ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai ${libraryDir}

echo $(date +"[%b %d %H:%M:%S]")

# Run ATACSeqQC R script
echo "Running ATACseqQC R script"

Rscript ${scriptDir}/ATACseqQC.R ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam ${outDir}/

echo $(date +"[%b %d %H:%M:%S]")

# Remove the bam and bai files from the ATACseqQC directory
echo "Removing bam and bai files from ATACseqQC directory"

rm ${libraryDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam
rm ${libraryDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

echo $(date +"[%b %d %H:%M:%S]")