#!/bin/bash

## Script for generating normalized bedgraphs
## Date: 19 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## bamDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## peakDir=/Users/coke6162/project20xx_cutnrun/macs2 \
## outDir=/Users/coke6162/project20xx_cutnrun/macs2/FRIP \
## expName=bosTau9.main \
## spikeName=Ecoli \
## sbatch --array=0-4 calculateFRIP.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=16G

# Job name and output
#SBATCH -J FRIP
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load modules
module load bedtools
module load samtools

# Define query files
queries=($(ls ${bamDir}/*.sorted.bam | xargs -n 1 basename | sed 's/_'${expName}'.sorted.bam//g' | sed 's/_'${spikeName}'.sorted.bam//g' | uniq))

# Set a variable equal to the total number of reads
pwd; hostname; date

echo "Bedtools version: "$(bedtools --version)
echo "Samtools version: "$(samtools --version)

echo "Processing file: "${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Setting variable equal to total number of reads..."

totalReads=`samtools view -c ${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${expName}.sorted.bam`

echo "Total number of reads: "$(echo $totalReads)

# Set a variable equal to the number of reads in peaks
echo $(date +"[%b %d %H:%M:%S] Setting a vairable equal to the number of reads in peaks...")

peakReads=`bedtools sort -i ${peakDir}/${queries[$SLURM_ARRAY_TASK_ID]}_peaks.narrowPeak \
| bedtools merge -i - \
| bedtools intersect -u -nonamecheck -a ${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${expName}.sorted.bam -b - -ubam \
| samtools view -c`

echo "Number of reads in peaks: "$(echo $peakReads)

# Calculate FRIP
echo $(date +"[%b %d %H:%M:%S] Calculating FRIP...")

FRIP=`awk "BEGIN {print "${peakReads}"/"${totalReads}"}"`

echo "${queries[$SLURM_ARRAY_TASK_ID]}" > ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Total number of mapped reads: "$(echo $totalReads) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Number of mapped reads in peaks: "$(echo $peakReads) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Fraction of mapped reads in peaks (FRIP): "$(echo $FRIP) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt

echo $(date +"[%b %d %H:%M:%S] Done")