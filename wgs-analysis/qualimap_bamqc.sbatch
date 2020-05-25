#!/bin/bash

## Script for running qualimap bamqc
## Date: 11 February 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=20200204_wgs_bosTau9_USMARC_BL3_hybrid \
## sbatch --array 0-6 qualimap_bamqc.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=30GB

# Job name and output
#SBATCH -J qualimap_bamqc
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Define path variables
inDir=${parentDir}/${projectName}/melt
outDir=${parentDir}/${projectName}/coverage/qualimap
tmpDir=${parentDir}/${projectName}/java_tmp

# Define query files
queries=($(ls ${inDir}/*_merged_marked_realigned_recalibrated.sorted.bam | xargs -n 1 basename))

# Make output directory
pwd; hostname; date

echo "Making output directory..."

mkdir -p ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged_marked_realigned_recalibrated.sorted.bam}

# Run qualimap
echo $(date +"[%b %d %H:%M:%S] Run qualimap bamqc...")

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

/Shares/CL_Shared/programs/qualimap_v2.2.1/qualimap \
bamqc \
-bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-outdir ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged_marked_realigned_recalibrated.sorted.bam} \
-outfile ${queries[$SLURM_ARRAY_TASK_ID]%.bam}.pdf \
-outformat PDF \
-nt ${numThreads} \
--java-mem-size=30G \
-c \
-nw 400 \
-hm 3

# Write file containing coverage info (for easy reference)
echo $(date +"[%b %d %H:%M:%S] Extracting coverage info...")

grep "mean coverageData" \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged_marked_realigned_recalibrated.sorted.bam}/genome_results.txt \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged_marked_realigned_recalibrated.sorted.bam}/mean_coverage.txt

echo $(date +"[%b %d %H:%M:%S] Done")