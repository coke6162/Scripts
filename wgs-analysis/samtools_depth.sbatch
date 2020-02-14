#!/bin/bash

## Script for running samtools depth
## Date: 11 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=20200204_wgs_bosTau9_USMARC_BL3_hybrid \
## sbatch --array 0-6 samtools_depth.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=30GB

# Job name and output
#SBATCH -J depth
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=1

# Define path variables
inDir=${parentDir}/${projectName}/melt
outDir=${parentDir}/${projectName}/coverage/samtools
tmpDir=${parentDir}/${projectName}/java_tmp

# Define query files
queries=($(ls ${inDir}/*_merged_marked_realigned_recalibrated.sorted.bam | xargs -n 1 basename))

# Load modules
module load samtools

# Calculate genome size
pwd; hostname; date

echo "samtools version: "$(samtools --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Calculating genome size..."

samtools view \
-H ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| grep -P '^@SQ' \
| cut -f 3 -d ':' \
| awk '{sum+=$1} END {print sum}' \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_genomeSize.txt

genomeSize=$(cat ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_genomeSize.txt)

# Run samtools
echo $(date +"[%b %d %H:%M:%S] Running samtools depth...")

samtools depth \
-a \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_samtools_depth.txt

# Calculate genome coverage
echo $(date +"[%b %d %H:%M:%S] Calculating genome depth...")

awk -v var="$genomeSize" '{sum+=$3} END { print "Average = ",sum/var}' \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_samtools_depth.txt \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_samtools_depth_genomeCov.txt

echo $(date +"[%b %d %H:%M:%S] Done")