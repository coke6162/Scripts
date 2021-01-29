#!/bin/bash

## Script for running tss enrichment R script
## Date: 29 January 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## scriptDir=/Users/coke6162/repos/scripts/atacseq-analysis \
## inDir=/scratch/Users/coke6162/project20xx_atacseq/bowtie2/tss_enrich \
## outDir=/scratch/Users/coke6162/project20xx_atacseq/bowtie2/tss_enrich \
## sbatch --array=0-3 tss_enrichment_score.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=1G

# Job name and output
#SBATCH -J tss_enrich
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*matrix | xargs -n 1 basename))

# Load modules
module load R/3.5.1

# Specify R library paths
export R_LIBS='/Users/coke6162/R/x86_64-redhat-linux-gnu-library/3.5/'

# Run R script with two paremeters, a path and the input matrix.
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

Rscript ${scriptDir}/tss_enrichment_score.R ${inDir}/ ${queries[$SLURM_ARRAY_TASK_ID]} ${outDir}/

echo $(date +"[%b %d %H:%M:%S] Done")