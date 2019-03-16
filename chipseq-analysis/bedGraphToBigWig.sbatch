#!/bin/bash

## Script for converting bedgraphs to bigwigs
## Date: 16 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/macs2 \
## outDir=/Users/coke6162/project20xx_chipseq/bigwig \
## sbatch --array 0-42 bedGraphToBigWig.sbatch 

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bedGraphToBigWig
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define key variables 
chromSizesFile=/scratch/Shares/chuong/genomes/mm10/mm10.chrom.sizes

# Other options include:
# /scratch/Shares/chuong/genomes/hg38/hg38.chrom.sizes 
# /scratch/Shares/chuong/genomes/hg19/hg19.chrom.sizes 
# /scratch/Shares/chuong/genomes/bosTau8/bosTau8.chrom.sizes

# Define query files
queries=($(ls ${inDir}/*treat*.bdg | xargs -n 1 basename))

# Convert bedgraphs to bigwigs
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bedgraph...")

sort -k1,1 -k2,2n ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} > ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Converting sorted bedgraph to bigwig...")
echo "bedGraphToBigWig v 4"

bedGraphToBigWig ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg ${chromSizesFile} ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.bw

echo $(date +"[%b %d %H:%M:%S] Removing sorted bedgraph...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Done!")