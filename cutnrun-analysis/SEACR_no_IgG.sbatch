#!/bin/bash

## Script for running SEACR
## Date: September 13 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## outDir=/Users/coke6162/project20xx_cutnrun/seacr/relaxed \
## expName=bosTau9.main \
## spikeName=Ecoli \
## threshold=1 \
## mode=relaxed \
## sbatch --array=0-4 SEACR_no_IgG.sbatch

## genome sizes not hardcoded: bosTau8 (2.66e9), canFam3 (2.33e09)

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J macs2
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Module load macs2
module load python/2.7.14/MACS/2.1.1

# Define query files
queries=($(ls ${inDir}/*.sam | xargs -n 1 basename | sed 's/_'${expName}'.sam//g' | sed 's/_'${spikeName}'.sam//g' | uniq))

# run the thing
pwd; hostname; date

echo "Calling peaks with SEACR..."

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo "Numeric threshold: "${threshold}

bash /Users/CL_Shared/programs/SEACR/SEACR_1.1.sh \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${expName}_d1000_fragments.bdg \
${threshold} \
non \
${mode} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${threshold}

echo $(date +"[%b %d %H:%M:%S] Done")