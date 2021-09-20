#!/bin/bash

## Script for converting bed files to fasta using bedtools getfasta
## Date: 15 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/input \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/input/fa \
## genomeFaFile=/Shares/CL_Shared/db/genomes/bosTau9/fa/bosTau9.fa \
## sbatch --array=0-0 bedtools_getfasta.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1G

# Job name and output
#SBATCH -J getfasta
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.bed | xargs -n 1 basename | sed 's/.bed//g' | uniq))

# load modules
module load bedtools

# make output directory
pwd; hostname; date

echo "Making output directory..."

mkdir -p ${outDir}

# run the thing
echo "Starting bedtools getfasta..."
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed

bedtools getfasta \
-fi ${genomeFaFile} \
-bed ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed \
-fo ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa

echo $(date +"[%b %d %H:%M:%S] Done")