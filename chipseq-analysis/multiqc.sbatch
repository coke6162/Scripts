#!/bin/bash

## Script for running multiqc
## Date: 31 Jan 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/fastqc \
## outDir=/Users/coke6162/project20xx_chipseq/fastqc/multiqc \
## sbatch multiqc.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J multiqc
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load FastQC
module load multiqc

# Define constant variables
numThreads=8

# Run MultiQC
pwd; hostname; date

echo "Running MultiQC..."
echo "MultiQC version: "(multiqc --version)

multiqc ${inDir}/*.zip --outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")

