#!/bin/bash

## Script for generating un-normalized bedgraphs
## Date: 13 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## outDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## expName=bosTau9.main \
## spikeName=Ecoli \
## genomeChrFile=/Users/CL_Shared/db/genomes/bosTau9/fa/bosTau9.main.chrom.sizes \
## sbatch --array=0-4 bamToBedGraph.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J bamToBdg
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Module load samtools and bedtools
module load samtools
module load bedtools

# Define query files
queries=($(ls ${inDir}/*.sam | xargs -n 1 basename | sed 's/_'${expName}'.sam//g' | sed 's/_'${spikeName}'.sam//g' | uniq))

# Generate un-normalized bedgraphs
pwd; hostname; date

echo "Generate un-normalized bedgraphs..."
echo "Samtools version: "$(samtools --version)
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

samtools sort \
-@ 8 \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${expName}.sorted.bam \
| bedtools genomecov \
-bg \
-ibam - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${expName}_unNorm.bdg

echo $(date +"[%b %d %H:%M:%S] Done")