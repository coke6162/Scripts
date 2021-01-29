#!/bin/bash

## Script for running deeptools to shift fragments +4 and -5 for positive and negative strand, respectively
## Date: 27 January 2021

## Example usage:
## cd /Users/coke6162/repos/project_files/yyyymmdd_ATAC-seq_projectName
## inDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/bowtie2/dedup \
## outDir=/scratch/Users/coke6162/yyyymmdd_ATAC-seq_projectName/bowtie2/dedup/shifted \
## sbatch --array 0-17 deeptools_alignmentSieve.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J alignmentSieve
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# load modules
module load singularity
module load samtools

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# Define query files
queries=($(ls ${inDir}/*.filtered.dedup.sorted.bam | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate enrichment matrix
numCPU=8

singularity exec --bind /scratch ${deeptools} \
alignmentSieve \
--bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--outFile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.filtered.dedup.sorted.bam}.filtered.dedup.shifted.sorted.bam.tmp \
--ATACshift \
--numberOfProcessors ${numCPU}

# Sort shifted bams
echo $(date +"[%b %d %H:%M:%S] Sorting shifted bams...")

samtools sort \
-@ ${numCPU} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.filtered.dedup.sorted.bam}.filtered.dedup.shifted.sorted.bam.tmp \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.filtered.dedup.sorted.bam}.filtered.dedup.shifted.sorted.bam

# Index sorted bam files
echo "Indexing sorted bams..."

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.filtered.dedup.sorted.bam}.filtered.dedup.shifted.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.filtered.dedup.sorted.bam}.filtered.dedup.shifted.sorted.bam.bai

# Remove intermediate bam files
echo "Removing intermediate bam files..."

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.filtered.dedup.sorted.bam}.filtered.dedup.shifted.sorted.bam.tmp

echo $(date +"[%b %d %H:%M:%S] Done")