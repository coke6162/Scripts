#!/bin/bash

## Script to download a bunch of sras
## Date: 16 July 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis/
## inDir=/Users/coke6162 \
## projectName=project20xx_rnaseq \
## sbatch --array 0-40 sraDownload.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=6:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J sraDownload
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -o /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load sra

# set query files
queries=($(cat ${inDir}/repos/project_files/${projectName}/sraDownload.txt | awk '(NR>1){print $1}'))

# run the thing
pwd; hostname; date

echo "Working directory: " ${inDir}/${projectName}/fastq
echo "fastq-dump version: "$(fastq-dump --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

cd ${inDir}/${projectName}/fastq

fastq-dump --split-files --gzip ${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done")