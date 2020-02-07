#!/bin/bash

## Script to download a bunch of sras
## Date: 3 February 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## sbatch --array 0-35 sraDownload.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J sraDownload
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load sra

# Define path variables
inDir=${parentDir}/${projectName}/fastq
refDir=${parentDir}/repos/project_files/${projectName}

# Define query files
# This works with a txt file where the first column (starting with row #2) contains SRR IDs
# If your SRR IDs are in a column other than the first, modify the 'print $1' argument
# If you do not have a header, remove the 'NR>1' argument
queries=($(cat ${refDir}/sraDownload.txt | awk '(NR>1){print $1}'))

# Run fastq-dump
pwd; hostname; date

echo "Working directory: " ${inDir}
echo "fastq-dump version: "$(fastq-dump --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

cd ${inDir}

fastq-dump --split-files --gzip ${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done")