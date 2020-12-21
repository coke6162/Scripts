#!/bin/bash

## Script to copy data from /scratch
## Date: 21 December 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## inDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_RNA-seq_projectName \
## destinationDir=/Users/coke6162 \
## sbatch copy_from_scratch.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=2:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J cp
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# run the thing
pwd; hostname; date

echo "Working dir: "${inDir}/${projectName}
echo "Copying to: "${destinationDir}

# make necessary directories
echo $(date +"[%b %d %H:%M:%S] Making destination directories...")

mkdir -p ${destinationDir}/${projectName}/fastq/fastqc/multiqc
mkdir -p ${destinationDir}/${projectName}/hisat2
mkdir -p ${destinationDir}/${projectName}/featureCounts
mkdir -p ${destinationDir}/${projectName}/bigwigs
mkdir -p ${destinationDir}/${projectName}/DESeq2

# go to the designated work dir
echo $(date +"[%b %d %H:%M:%S] Moving to working directory...")

cd ${inDir}/${projectName}

# make necessary subdirs
cp -n -r fastq/ ${destinationDir}/${projectName}/
cp -n -r hisat2/ ${destinationDir}/${projectName}/
cp -n -r featureCounts/ ${destinationDir}/${projectName}/
cp -n -r bigwigs/ ${destinationDir}/${projectName}/
cp -n -r DESeq2/ ${destinationDir}/${projectName}/

echo $(date +"[%b %d %H:%M:%S] Done")