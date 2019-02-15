#!/bin/bash

## Script for running trimmomatic
## Date: 31 Jan 2019 

## Example usage:
## inDir=/Users/%u/project20xx_chipseq/fastq \
## outDir=/Users/%u/project20xx_chipseq/trimmomatic \
## logOut=/Users/%u/project20xx_chipseq/trimmomatic/log \
## adapterFile=/opt/trimmomatic/0.36/adapters/TruSeq3-PE.fa \
## sbatch --array 0-42 trimmomatic_PE.sbatch

## adapterFile options: NexteraPE-PE.fa, TruSeq2-PE.fa, TruSeq3-PE-2.fa, TruSeq3-PE.fa

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=8:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J trimmomatic
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@colorado.edu

# Set constant variables
trimExe=/opt/trimmomatic/0.36/trimmomatic-0.36.jar
numThreads=8

# Define query files
# Removes R1/R2 to generate a unique identifier for each pair of files
# Modify to fit input file name if necessary
queries=($(ls ${inDir}/*fastq.gz | xargs -n 1 basename | sed --expression='s/_R1.fastq.gz//g' \
	| sed --expression='s/_R2.fastq.gz//g' | uniq))

# Run Trimmomatic
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Trimmomatic version: "$(java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar --version)
echo "Starting trimmomatic..."

java -jar ${trimExe} PE -threads ${numThreads} \
-trimlog ${logOut}/${queries[$SLURM_ARRAY_TASK_ID]}_trimlog.txt \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1.fastq.gz ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2.fastq.gz \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_unpaired.fastq.gz \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_unpaired.fastq.gz \
ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/${adapterFile}:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# 'PE' - paired-end input
# '-threads <int>' - number of threads to use
# '-trimlog <file.txt>' - creates a log of all read trimmings 
# 'ILLUMINACLIP <file.fa>' - cut adapter and other illumina-specific sequences from the read
# 'SLIDINGWINDOW <int:int>' - scans 5' end of reads and clips once average quality within the specified window falls below a threshold
# 'LEADING <int>' - cut bases off the start of a read, if below a threshold quality
# 'TRAILING <int>' - cut bases off the end of a read, if below a threshold quality
# 'MINLEN <int>' - retain only reads of some length (in bp)