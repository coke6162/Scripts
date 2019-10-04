#!/bin/bash

## Script for running salmon quant with PE reads
## Date: 21 April 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## indexDir=/Users/CL_Shared/db/genomes/mm10/index/salmon/gencode.vM18.transcripts_salmon_index \
## gtfFile=/Users/CL_Shared/db/genomes/mm10/annotations/gencode.vM18.annotation.gtf \
## inDir=/Users/coke6162/projoect20xx_rnaseq/trimmomatic \
## outDir=/Users/coke6162/project20xx_rnaseq/salmon \
## sbatch salmon_quant_PE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J salmon_quant_PE
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
queries=($(ls ${inDir}/*_paired.fastq.gz | xargs -n 1 basename | sed 's/_R1_paired.fastq.gz//g' \
	| sed 's/_R2_paired.fastq.gz//g' | uniq))

# Load hisat2 and samtools
module load salmon

# Run hisat2
pwd; hostname; date

echo "salmon version: "$(salmon --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Starting salmon pseudo-alignment..."

salmon quant -i ${indexDir} --libType A -p ${numThreads} --geneMap ${gtfFile} \
-1 <(gunzip -c ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz) \
-2 <(gunzip -c ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz) \
--validateMappings --rangeFactorizationBins 4 --gcBias \
-o ${outDir}/{queries[$SLURM_ARRAY_TASK_ID]}_transcripts_quant

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-i <indexDir>' - location of salmon index directory
# '-libType <arg>' - library type; specify '-libType A' for automatic detection or see manual for other options
# '-p <int>' - number of threads
# '--geneMap <gtfFile>' - gtf file containing mapping of transcripts to genes; allows for gene-level abundance estimates
# '-1 <file>' - file containing #1 mates
# '-2 <file>' - file containing #2 mates
# '--validateMappings' - validate mappings using alignment-based verification (highly recommended)
# '--rangeFactorizationBins <int>' - factorizes the likelihood used in quantification; common value to select is 4
# '--gcBias' - perform GC bias correction; typically recommended
# '-o <dir>' - ooutput directory