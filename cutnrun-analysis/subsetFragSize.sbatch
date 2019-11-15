#!/bin/bash

## Script for running deeptools 
## Date: 10 October 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/20190910_CUTnRUN_MDBK/bwa/filtered \
## outDir=/Users/coke6162/20190910_CUTnRUN_MDBK/bwa/filtered \
## sbatch --array 0-12 subsetFragSize.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J subsetFragSize
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.bam | xargs -n 1 basename | sed --expression='s/_filtered.sorted.bam//g' | uniq))

# load modules
module load singularity samtools

# define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
numThreads=8

# Use "alignmentSieve" to filter reads by size (<120 or >120bp)
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Starting deeptools alignmentSieve...")

singularity exec --bind /Shares/CL_Shared ${deeptools} \
alignmentSieve \
--bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--maxFragmentLength 120 \
--outFile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.bam.tmp \
--filteredOutReads ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.bam.tmp \
--numberOfProcessors ${numThreads} \

# Sort and index output bams
echo $(date +"[%b %d %H:%M:%S] Sorting bams (less than 120)...")

samtools sort \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.bam.tmp \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.sorted.bam

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.bam.tmp

echo $(date +"[%b %d %H:%M:%S] Sorting bams (more than 120)...")

samtools sort \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.bam.tmp \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.sorted.bam

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.bam.tmp

echo $(date +"[%b %d %H:%M:%S] Indexing bams (less than 120)...")

samtools index \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.sorted.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.sorted.bam.bai

echo $(date +"[%b %d %H:%M:%S] Indexing bams (more than 120)...")

samtools index \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.sorted.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.sorted.bam.bai

# Generate flagstat files from bam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from bams (all fragments)...")

samtools flagstat \
-@ ${numThreads} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_flagstat.txt

echo $(date +"[%b %d %H:%M:%S] Making flagstat file from bams (less than 120)...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered.sorted.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan120_filtered_flagstat.txt

echo $(date +"[%b %d %H:%M:%S] Making flagstat file from bams (more than 120)...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered.sorted.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_moreThan120_filtered_flagstat.txt

echo $(date +"[%b %d %H:%M:%S] Done")