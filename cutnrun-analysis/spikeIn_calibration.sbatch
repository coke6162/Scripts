#!/bin/bash

## Script for generating normalized bedgraphs
## Date: 11 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## outDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## expName=bosTau9.main \
## spikeName=Ecoli \
## genomeChrFile=/Users/CL_Shared/db/genomes/bosTau9/fa/bosTau9.main.chrom.sizes \
## sbatch --array=0-4 spikeIn_calibration.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J spikeIn_cal
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Module load bowtie2 and samtools
module load bedtools

# Define query files
queries=($(ls ${inDir}/*.sam | xargs -n 1 basename | sed 's/_'${expName}'.sam//g' | sed 's/_'${spikeName}'.sam//g' | uniq))

# Determine normalization factor (number of reads mapping to E. coli per 10000 reads) used to generate normalized bedgraph
pwd; hostname; date

echo "Aligning reads to the genome & outputting sam files"
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Number of reads mapping to bosTau9: "$(wc -l ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_bosTau9.main_d1000_fragments.sorted.bed)
echo "Number of reads mapping to E. coli: "$(wc -l ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_Ecoli_d1000_fragments.sorted.bed)

spikeCount=`wc -l ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_Ecoli_d1000_fragments.sorted.bed | awk '{print $1}'`
normFactor=`echo "10000 / ${spikeCount}" | bc -l`
echo "normFactor =" ${normFactor}

# Generate normalized bedgraphs
echo "Generating normalized bedgraphs..."
sort -k1,1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_bosTau9.main_d1000_fragments.sorted.bed \
| bedtools genomecov \
-bg \
-scale ${normFactor} \
-i - \
-g ${genomeChrFile} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_bosTau9.main_d1000_fragments.bdg

echo $(date +"[%b %d %H:%M:%S] Done")