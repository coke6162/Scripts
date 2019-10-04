#!/bin/bash

## Script for running hisat2 index
## Date: 4 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## refFile=/Users/CL_Shared/db/genomes/bosTau8/fa/bosTau8.fa \
## outDir=/Users/CL_Shared/db/genomes/bosTau8/index \
## name=bosTau8 \
## sbatch hisat2_index.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=6:00:00
#SBATCH --mem=32G

# Job name and output
#SBATCH -J hisat2_index
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load module
module load hisat2

# Run hisat2-build
pwd; hostname; date

echo "hisat2-build version: "$(hisat2-build --version)
echo $(date +"[%b %d %H:%M:%S] Running hisat2-build...")

hisat2-build -f --threads ${numThreads} ${refFile} ${outDir}/${name}

echo $(date +"[%b %d %H:%M:%S] Done")