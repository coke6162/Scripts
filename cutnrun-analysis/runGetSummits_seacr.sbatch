#!/bin/bash

## Script for running getSummits_seacr.py as an array
## Date: 18 November 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/20191111_CUTnRUN_MDCK/seacr/stringent \
## outDir=/Users/coke6162/20191111_CUTnRUN_MDCK/seacr/stringent \
## mode=stringent \
## sbatch --array=0-38 runGetSummits_seacr.sbatch

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
queries=($(ls ${inDir}/*.${mode}.bed | xargs -n 1 basename | sed --expression='s/.'${mode}'.bed//g' | uniq))

# run the thing
pwd; hostname; date

echo "Starting getSummits_seacr.py..."

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed

/opt/python/3.6.3/bin/python3.6 getSummits_seacr.py \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.${mode}.bed \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.${mode}_summits.bed

echo $(date +"[%b %d %H:%M:%S] Done")