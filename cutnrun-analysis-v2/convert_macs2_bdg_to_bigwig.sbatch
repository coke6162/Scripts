#!/bin/bash

## Script for converting macs2 bedgraphs to bigwigs
## Date: 11 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/5_macs2_output_PE \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/7_bigwigs \
## chromSizes=/Shares/CL_Shared/db/genomes/bosTau9/fa/bosTau9.chrom.sizes \
## sbatch --array 0-0 convert_macs2_bdg_to_bigwig.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=2GB

# Job name and output
#SBATCH -J bedGraphToBigWig
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*_treat_pileup.bdg | xargs -n 1 basename | sed 's/_treat_pileup.bdg//g'))

# Sort bedgraph
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bedgraph...")

sort -k1,1 -k2,2n ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_treat_pileup.bdg \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_macs2.bdg.tmp

# Convert sorted bedgraph to bigwig
echo $(date +"[%b %d %H:%M:%S] Converting sorted bedgraph to bigwig...")

/Shares/CL_Shared/programs/ucsc/bedGraphToBigWig \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_macs2.bdg.tmp \
${chromSizes} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_macs2.bw

# Remove sorted bedgraph
echo $(date +"[%b %d %H:%M:%S] Removing sorted bedgraph...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_macs2.bdg.tmp

echo $(date +"[%b %d %H:%M:%S] Done")