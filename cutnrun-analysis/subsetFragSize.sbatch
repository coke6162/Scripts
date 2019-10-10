#!/bin/bash

## Script for running deeptools 
## Date: 10 October 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/20190910_CUTnRUN_MDBK/bwa/filtered \
## outDir=/Users/coke6162/20190910_CUTnRUN_MDBK/bwa/filtered/subsetFragSize \
## sbatch --array 0-12 subsetFragSize.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1:00:00
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
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename | uniq))

# load modules
module load singularity

# define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
numThreads=8

# Use "alignmentSieve" to filter reads by size (<120 or >120bp)
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Starting deeptools alignmentSieve...")

singularity exec --bind /Shares/CL_Shared ${deeptools} \
alignmentSieve \
--bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--maxFragmentLength 120 \
--outFile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_small.bam \
--filterMetrics ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_small_metrics.txt \
--filteredOutReads ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_large.bam
--numberOfProcessors ${numThreads} \

echo $(date +"[%b %d %H:%M:%S] Done")