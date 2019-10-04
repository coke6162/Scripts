#!/bin/bash

## Script for running sambamba markdup
## Date: 14 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/bwa \
## outDir=/Users/coke6162/project20xx_chipseq/bwa \
## sbatch --array 0-42 sambamba_markdup.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J sambamba_markdup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load samtools
module load samtools

# Define query files
queries=($(ls ${inDir}/*.sorted.uniq.bam | xargs -n 1 basename))

# Remove duplicates with sambamba
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running sambamba markdup...")

sambamba markdup --remove-duplicates --nthreads=4 --overflow-list-size=600000 --tmpdir=${inDir}/tmp --show-progress \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.uniq.bam}_dedup.bam

echo $(date +"[%b %d %H:%M:%S] Done")

# Run sambamba flagstat on the input and output to get a sense of removed reads
pwd; hostname; date

echo "Processing input file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Running sambamba flagstat on input bam"

sambamba flagstat ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} > ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.uniq.bam}_flagstat.txt

echo "Processing dedup file: "${queries[$SLURM_ARRAY_TASK_ID]%.sorted.uniq.bam}_dedup.bam
echo "Running sambamba flagstat on dedup bam"

sambamba flagstat ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.uniq.bam}_dedup.bam > ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.uniq.bam}_dedup_flagstat.txt

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '--remove-duplicates' - remove duplicates instead of just marking them
# '-nthreads=<int>' - number of threads to use
# '--overflow-list-size=<int>' - size of the overflow list where reads, thrown from the hash table, get a second chance to meet their pairs ( ͡° ͜ʖ ͡°)
# '--tmpdir=<tmpDir>' - specify directory for temporary files