#!/bin/bash

## Script for merging macs2 narrowPeak files from SE and PE modes
## Date: 11 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir1=/scratch/Users/coke6162/yyyymmdd_projectName/4_macs2_output_SE/histone \
## inDir2=/scratch/Users/coke6162/yyyymmdd_projectName/5_macs2_output_PE/histone \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/6_macs2_output_merged/histone  \
## sbatch --array=0-0 merge_peak_files.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1G

# Job name and output
#SBATCH -J merge_peaks
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load bedtools

# Define query files
queries=($(ls ${inDir1}/*_peaks.narrowPeak | xargs -n 1 basename))

# Run the thing
pwd; hostname; date

echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir1}/${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Merging peak files...")

cat ${inDir1}/${queries[$SLURM_ARRAY_TASK_ID]} ${inDir2}/${queries[$SLURM_ARRAY_TASK_ID]} \
| bedtools sort -i - \
| bedtools merge -i - -c 5 -o max \
| awk '{print $1 "\t" $2 "\t" $3 "\t" "'${queries[$SLURM_ARRAY_TASK_ID]%_peaks.narrowPeak}_peak'"NR "\t" $4 "\t" "."}' \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_peaks.narrowPeak}_mergedpeaks.bed

echo $(date +"[%b %d %H:%M:%S] Done!")
