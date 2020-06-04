#!/bin/bash

## Script for download fastq files (via SRA identifier)
## Date: 1 July 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## parentDir=/Users/coke6162 \
## projectName=megquier2019_ATAC-seq_cf \
## sbatch sraDownlaod.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J sraDownload
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load sra

# set query files
queries=($(cat ${parentDir}/repos/project_files/${projectName}/sraDownload.txt \
	| awk '(NR>1){print $1}'))

# run the thing
pwd; hostname; date

echo "Working directory: "${parentDir}/repos/project_files/${projectName}/
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

cd ${parentDir}/${projectName}/fastq

fastq-dump --split-files --gzip ${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done")