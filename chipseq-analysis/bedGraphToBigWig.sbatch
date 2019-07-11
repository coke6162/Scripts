#!/bin/bash

## Script for converting bedgraphs to bigwigs
## Date: 11 July 2019 

## Example usage:
## cd /Users/coke6162/repos/scripts/chipseq-analysis
## inDir=/Users/coke6162/project20xx_chipseq/macs2 \
## outDir=/Users/coke6162/project20xx_chipseq/bigwigs \
## chromSizesFile=/Users/CL_Shared/db/genomes/mm10/fa/mm10.chrom.sizes \
## sbatch --array 0-42 bedGraphToBigWig.sbatch 

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
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
queries=($(ls ${inDir}/*treat*.bdg | xargs -n 1 basename))

# Convert bedgraphs to bigwigs
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bedgraph...")

sort -k1,1 -k2,2n \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Converting sorted bedgraph to bigwig...")
echo "bedGraphToBigWig v 4"

bedGraphToBigWig \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg \
${chromSizesFile} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.bw

echo $(date +"[%b %d %H:%M:%S] Removing sorted bedgraph...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Done")