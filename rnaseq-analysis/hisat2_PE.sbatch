#!/bin/bash

## Script for running hisat2
## Date: 30 Narch 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_rnaseq/bbduk \
## outDir=/Users/coke6162/project20xx_rnaseq/hisat2 \
## hisatIndexDir=/Users/CL_Shared/db/genomes/bosTau8/index/hisat2 \
## hisatIndex=bosTau8 \
## rnaStrandness=FR \
## sbatch --array 0-42 hisat2_PE.sbatch

## hisatIndexDir options: hg38, mm10, bosTau8
## hisatIndex: hg38.main, mm10, bosTau8
## rnaStrandness options: FR (forward) or RF (reverse stranded)

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J hisat2_PE
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls ${inDir}/*_bbduk.fastq.gz | xargs -n 1 basename | sed 's/_1_bbduk.fastq.gz//g' \
	| sed 's/_2_bbduk.fastq.gz//g' | uniq))

# Load hisat2 and samtools
module load hisat2 samtools

# Run hisat2
pwd; hostname; date

echo "HISAT2 version: "$(hisat2 --version)
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running hisat2...")

hisat2 \
-p ${numThreads} \
--rna-strandness ${rnaStrandness} \
--no-softclip \
-x ${hisatIndexDir}/${hisatIndex} \
-1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1_bbduk.fastq.gz \
-2 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2_bbduk.fastq.gz \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam

# Create .sorted.uniq.bai index files for upload to UCSC browser
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-p <int>' - number of threads
# '--rna-strandness <rf/fr>' - specify strand-specific information; library kit typically has this information
# '--no-softclip' - disallows soft clipping
# '-x <file.basename>' - hisat2 index basename; include everything up to .main.1.ht2
# '-1 <file_R1_bbduk.fastq.gz' - mate file 1
# '-2 <file_R2_bbduk.fastq.gz' - mate file 2
# '-q <int>' - filter out mapped reads by quality threshold
# '-Sb' - convert sam to bam
# '-@' - number of threads