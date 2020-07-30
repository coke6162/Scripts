#!/bin/bash

## Script to set up directory layout 
## Date: 4 May 2019

## Example usage:
## inDir=/Users/coke6162 \
## projectName=project20xx_rnaseq \
## sbatch setup.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J setup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# run the thing
pwd; hostname; date

echo "Working dir: "${inDir}
echo "Project: "${projectName}
echo $(date +"[%b %d %H:%M:%S] Making project directories...")

# go to the designated work dir
cd ${inDir}

# make a new project dir
mkdir -p ${projectName}

# make necessary subdirs
mkdir -p ${projectName}/fastq
mkdir -p ${projectName}/fastq/fastqc
mkdir -p ${projectName}/fastq/fastqc/multiqc
mkdir -p ${projectName}/bbduk
mkdir -p ${projectName}/bbduk/log
mkdir -p ${projectName}/bbduk/fastqc
mkdir -p ${projectName}/bbduk/fastqc/multiqc
mkdir -p ${projectName}/hisat2
mkdir -p ${projectName}/featureCounts
mkdir -p ${projectName}/bigwigs
mkdir -p ${projectName}/DESeq2
mkdir -p ${projectName}/DESeq2/png
mkdir -p ${projectName}/DESeq2/bed

echo $(date +"[%b %d %H:%M:%S] Done")