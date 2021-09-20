#!/bin/bash

## Script for calculating FRIP score 
## Date: 11 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## bamDir=/scratch/Users/coke6162/yyyymmdd_projectName/3_bams/histone \
## peakDir=/scratch/Users/coke6162/yyyymmdd_projectName/6_macs2_output_merged/histone \
## reportsDir=/scratch/Users/coke6162/yyyymmdd_projectName/reports/frip_scores \
## sbatch --array=0-0 calculate_frip_score.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=2:00:00
#SBATCH --mem=4G

# Job name and output
#SBATCH -J frip_score
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=1

# Load modules
module load bedtools samtools

# Define query files
queries=($(ls ${bamDir}/*_filtered.sorted.bam | xargs -n 1 basename | sed 's/_filtered.sorted.bam//g'))

# Run the thing
pwd; hostname; date

echo "Bedtools version: "$(bedtools --version)
echo "Samtools version: "$(samtools --version)

echo "Processing file: "${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Calculating total reads and peak reads...")

# Set a variable equal to the total number of reads (from bam file)
totalReads=`samtools view -c ${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam`

echo "Total number of reads: "$(echo $totalReads)

# Set a variable equal to the number of reads in peaks (from macs2 narrowpeak file)
peakReads=`bedtools sort -i ${peakDir}/${queries[$SLURM_ARRAY_TASK_ID]}_mergedpeaks.bed \
| bedtools merge -i - \
| bedtools intersect -u -nonamecheck -a ${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam -b - -ubam \
| samtools view -c`

echo "Number of reads in peaks: "$(echo $peakReads)

# Calculate FRiP score
echo $(date +"[%b %d %H:%M:%S] Calculating FRIP...")

FRIP=`awk "BEGIN {print "${peakReads}"/"${totalReads}"}"`

echo "${queries[$SLURM_ARRAY_TASK_ID]}" > ${reportsDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Total number of mapped reads: "$(echo $totalReads) >> ${reportsDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Number of mapped reads in peaks: "$(echo $peakReads) >> ${reportsDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Fraction of mapped reads in peaks (FRiP): "$(echo $FRIP) >> ${reportsDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt

echo $(date +"[%b %d %H:%M:%S] Done")