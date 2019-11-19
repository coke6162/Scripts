#!/bin/bash

## Script for running bedtools slop
## Date: 18 November 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/20191111_CUTnRUN_MDCK/seacr/stringent \
## outDir=/Users/coke6162/20191111_CUTnRUN_MDCK/seacr/stringent/slopSummits \
## genomeChrFile=/Shares/CL_Shared/db/genomes/canFam3/fa/canFam3.main.chrom.sizes \
## slopValue=75 \
## sbatch --array=0-38 bedtools_slop.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J runGetSummits_seacr.sbatch
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*_summits.bed | xargs -n 1 basename | sed --expression='s/_summits.bed//g' | uniq))

# load modules
module load bedtools

# run the thing
pwd; hostname; date

echo "Starting bedtools slop..."
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_summits.bed

bedtools slop \
-i ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_summits.bed \
-g ${genomeChrFile} \
-b ${slopValue} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_slopSummits_b${slopValue}.bed

echo $(date +"[%b %d %H:%M:%S] Done")