#!/bin/bash

## Script for running bbduk
## Date: 16 July 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## inDir=/Users/coke6162/project20xx_rnaseq/fastq \
## outDir=/Users/coke6162/project20xx_rnaseq/bbduk \
## adapterFile=/opt/bbmap/38.05/resources/adapters.fa \
## sbatch --array 0-42 bbduk_SE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=32G

# Job name and output
#SBATCH -J bbduk_SE
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load module
module load bbmap

# Define query files
queries=($(ls ${inDir}/*fastq.gz | xargs -n 1 basename))

# Run bbduk
pwd; hostname; date

echo "bbduk version: "$(bbduk.sh --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running bbduk...")

bbduk.sh -Xmx4g in=${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
out=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
ref=${adapterFile} \
ktrim=r k=34 mink=11 hdist=1 tpe tbo \
qtrim=r trimq=10 \
t=${numThreads}

echo $(date +"[%b %d %H:%M:%S] Done")