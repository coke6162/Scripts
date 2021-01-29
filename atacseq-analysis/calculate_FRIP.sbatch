#!/bin/bash

## Script for generating normalized bedgraphs
## Date: 27 January 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## refDir=/Users/coke6162/repos/project_files/yyyymmdd_ATAC-seq_projectName \
## bamDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/bowtie2/dedup\
## peakDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/macs2 \
## outDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/macs2/FRIP \
## sbatch --array=0-4 calculate_FRIP.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J FRIP
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load bedtools
module load samtools

# Define query samples
queries=($(cat ${refDir}/samples.txt))

# Specify names of bam and peak files
pwd; hostname; date

echo "Bedtools version: "$(bedtools --version)
echo "Samtools version: "$(samtools --version)

echo "Processing sample: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Setting variable to bam and peak files...")

bam_file=$(ls ${bamDir} | grep "^${queries[$SLURM_ARRAY_TASK_ID]}.*bam$" - | sed 's|^|'${bamDir}'/|')
peak_file=$(ls ${peakDir} | grep "^${queries[$SLURM_ARRAY_TASK_ID]}.*_peaks.narrowPeak$" | sed 's|^|'${peakDir}'/|')

# Set a variable equal to the total number of reads
echo $(date +"[%b %d %H:%M:%S] Setting a vairable equal to the number of reads in peaks...")

totalReads=`samtools view -c ${bam_file}`

echo "Total number of reads: "$(echo $totalReads)

# Calculating number of reads in peaks
echo $(date +"[%b %d %H:%M:%S] Calculating number of reads in peaks...")

peakReads=`bedtools sort -i ${peak_file} \
	| bedtools merge -i - \
	| bedtools intersect -u -nonamecheck -a ${bam_file} -b - -ubam \
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
