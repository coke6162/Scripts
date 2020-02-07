#!/bin/bash

## Script for converting bedgraphs to bigwigs
## Date: 11 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/project20xx_cutnrun/macs2 \
## outDir=/Users/coke6162/project20xx_cutnrun/bigwig \
## chromSizes=/Shares/CL_Shared/db/genomes/bosTau8/fa/bosTau9.chrom.sizes \
## sbatch --array 0-3 bedGraphToBigWig.sbatch 

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bedGraphToBigWig
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.bdg | xargs -n 1 basename))

# Convert bedgraphs to bigwigs
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bedgraph...")

sort \
-k1,1 -k2,2n \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Converting sorted bedgraph to bigwig...")
echo "bedGraphToBigWig v 4"

/Shares/CL_Shared/programs/ucsc/bedGraphToBigWig \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg \
${chromSizes} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.bw

echo $(date +"[%b %d %H:%M:%S] Removing sorted bedgraph...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Done")