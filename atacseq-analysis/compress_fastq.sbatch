#!/bin/bash

## Script to download a bunch of sras
## Date: 30 January 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## workingDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/fastq \
## referenceFile=/Users/coke6162/repos/project_files/yyyymmdd_ATAC-seq_projectName/SRR.txt \
## sbatch --array 0-17 compress_fastq.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J gzip
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
# This works with a txt file where the first column (starting with row #2) contains SRR IDs
# If your SRR IDs are in a column other than the first, modify the 'print $1' argument
# If you do not have a header, remove the 'NR>1' argument
queries=($(cat ${referenceFile} | awk '(NR>1){print $1}'))

# Compress fastq files
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Compressing fastq files...")

echo "Working directory: "${workingDir}
echo "Processing sample: "${queries[$SLURM_ARRAY_TASK_ID]}

gzip ${workingDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1.fastq
gzip ${workingDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2.fastq

echo $(date +"[%b %d %H:%M:%S] Done")