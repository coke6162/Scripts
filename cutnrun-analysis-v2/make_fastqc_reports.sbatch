#!/bin/bash

## Script for running fastqc
## Date: 10 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/1_raw_fastq \
## reportsDir=/scrath/Users/coke6162/yyyymmdd_projectName/reports/raw_fastqc \
## sbatch --array 0-7 make_fastqc_reports.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J fastqc
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Load FastQC
module load fastqc

# Define constant variables
numThreads=1

# Define query files
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename))

# Run FastQC
pwd; hostname; date

echo "FastQC version: "$(fastqc --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running fastqc...")

fastqc -o ${outDir} -f fastq -t ${numThreads} ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-f <file_type>' - input file type; options: fastq, bam, sam, bam_mapped, and sam_mapped
# '-t <int>' - number of threads