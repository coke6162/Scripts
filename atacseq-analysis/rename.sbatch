#!/bin/bash

## Script to rename fastq files
## Date: 26 January 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## workingDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/fastq \
## referenceFile=/Users/coke6162/repos/project_files/yyyymmdd_ATAC-seq_projectName/SRR.txt \
## sbatch --array 0-17 rename.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J rename
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# define query
# one job will spawn per SRR ID in reference file
queries=($(cat ${referenceFile} | awk '(NR>1){print $1}'))

# set variable for renaming
# this assumes column 4 is set to unique sample name
newName=$(grep ${queries[$SLURM_ARRAY_TASK_ID]} ${referenceFile} | awk -F " " '{print $2}')

# rename fastq files
pwd; hostname; date

echo "Working dir: "${workingDir}
echo "Renaming SRR: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Renaming fastq files...")

mv ${workingDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz ${workingDir}/${newName}_1.fastq.gz
mv ${workingDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2.fastq.gz ${workingDir}/${newName}_2.fastq.gz

echo $(date +"[%b %d %H:%M:%S] Done")