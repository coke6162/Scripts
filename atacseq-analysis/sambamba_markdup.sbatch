#!/bin/bash

## Script for running sambamba markdup
## Date: 1 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## inDir=/Users/coke6162/project20xx_atacseq/bowtie2 \
## outDir=/Users/coke6162/project20xx_atacseq/bowtie2/dedup \
## refFile=/Users/CL_Shared/db/genomes/canFam3/fa/canFam3.fa \
## sbatch --array=0-3 sambamba_markdup.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=4G

# Job name and output
#SBATCH -J sambamba_markdup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load modules
module load samtools
module load sambamba

# Define query files
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename))

# Make temporary directories
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${inDir}/tmp_sambamba
mkdir -p ${inDir}/tmp_bam/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}

# Filter out reads mapping to mitochondrial genes
echo $(date +"[%b %d %H:%M:%S] Removing fragments mapping to mitochondrial DNA")

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

samtools view \
-@ ${numThreads} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| egrep -v chrM \
| samtools view -bT ${refFile} -o ${inDir}/tmp_bam/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam

# Sort filtered bam files
echo $(date +"[%b %d %H:%M:%S] Sorting filtered bam files with sambamba")

sambamba sort \
--natural-sort \
--nthreads=${numThreads} \
--memory-limit=2G \
--tmpdir=${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam} \
--out=${inDir}/tmp_bam/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam \
${inDir}/tmp_bam/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam

# Remove duplicates with sambamba
echo $(date +"[%b %d %H:%M:%S] Removing duplicates with sambamba")

sambamba markdup \
--remove-duplicates \
--nthreads=${numThreads} \
--overflow-list-size=600000 \
--hash-table-size=4194304 \
--tmpdir=${inDir}/tmp \
${inDir}/tmp_bam/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.dedup.bam

# Run samtools flagstat 
echo $(date +"[%b %d %H:%M:%S] Running samtools flagstat on dedup'd bams")

samtools flagstat \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.dedup.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.dedup.flagstat.txt

# Sort dedup'd bams
echo $(date +"[%b %d %H:%M:%S] Sorting dedup'd bams")

samtools sort \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.dedup.bam \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.dedup.sorted.bam

# Index sorted dedup'd bams
echo $(date +"[%b %d %H:%M:%S] Indexing sorted dedup'd bams")

samtools index \
-@ {numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.dedup.sorted.bam

# Remove temporary files and directories
echo $(date +"[%b %d %H:%M:%S] Removing temporary files & directories")

rm -r ${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.filtered.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done")