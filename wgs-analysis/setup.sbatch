#!/bin/bash

## Script to set up directory structure
## Date: 3 February 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## inDir=/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## sbatch setup.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=0:05:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J fastqc
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set up directory structure
pwd; hostname; date

echo "Working dir: "${inDir}
echo "Project: "${projectName}
echo $(date +"[%b %d %H:%M:%S] Making project directories...")

# Go to the designated work dir
cd ${inDir}

# Make a new project dir
mkdir -p ${projectName}

# Make necessary subdirs
mkdir -p ${projectName}/fastq/fastqc/multiqc
mkdir -p ${projectName}/bbduk/fastqc/multiqc
mkdir -p ${projectName}/java_tmp
mkdir -p ${projectName}/bwa/flagstat
mkdir -p ${projectName}/bwa/picard/sorted
mkdir -p ${projectName}/bwa/picard/merged
mkdir -p ${projectName}/bwa/picard/marked/metrics
mkdir -p ${projectName}/gatk/realigned
mkdir -p ${projectName}/gatk/recalibrated
mkdir -p ${projectName}/melt
mkdir -p ${projectName}/coverage/samtools
mkdir -p ${projectName}/coverage/qualimap