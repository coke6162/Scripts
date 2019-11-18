#!/bin/bash

## Script for running SEACR
## Date: October 1 2019

## Example usage:
## cd /Users/coke6162/repos/project_files/20190910_CUTnRUN_MDBK
## inDir=/Users/coke6162/project20xx_cutnrun/bwa \
## outName=/Users/coke6162/project20xx_cutnrun/seacr/relaxed \
## mode=relaxed \
## bash runSEACR.sh

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J SEACR
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load bedtools R/3.5.1

# run the thing
pwd; hostname; date

echo "Calling peaks with SEACR..."
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir}/${inFile}

bash /Shares/CL_Shared/programs/SEACR_1.2/SEACR_1.2.sh \
${inDir}/${inFile} \
${inDir}/${controlFile} \
norm \
${mode} \
${outDir}/${outName}

echo $(date +"[%b %d %H:%M:%S] Done")