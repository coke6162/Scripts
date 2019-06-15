#!/bin/bash

## Script for running fimo
## Date: June 14 2019

## Example usage:
## cd /Users/coke6162/repos/project_files/20190613_ATAC-seq_BL3.1
## genomeName=bosTau9 \
## pvalue=1e-4 \
## motifFile=/Users/coke6162/repos/project_files/20190613_ATAC-seq_BL3.1/MA0137.3.meme \
## name=MA0137.3 \
## outDir=/Users/CL_Shared/db/genomes/bosTau9/fimo/MA0137.3 \
## sbatch fimo.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J fimo
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# load modules
module load meme
module load python/2.7.14/numpy/1.14.1
module load bedtools

# Define key variables
fastaFile=/Users/CL_Shared/db/genomes/${genomeName}/fa/${genomeName}.main.fa
chromSizes=/Users/CL_Shared/db/genomes/${genomeName}/fa/${genomeName}.main.chrom.sizes

# Run fimo
pwd; hostname; date

echo "fimo version: "$(fimo --version)
echo "bedtools version: "$(bedtools --version)

echo $(date +"[%b %d %H:%M:%S] Starting fimo...")

fimo \
--text \
--parse-genomic-coord \
--thresh ${pvalue} \
${motifFile} \
${fastaFile} \
> ${outDir}/pval_${pvalue}/${name}_${pvalue}.txt

# Convert output to bed & output to a new directory
echo $(date +"[%b %d %H:%M:%S] Convert fimo output to bed format...")

mkdir -p ${outDir}/pval_${pvalue}

/opt/python/2.7.14/bin/python2.7 /Users/CL_Shared/scripts/ed/fimo_analysis.py \
-i ${outDir}/pval_${pvalue}/${name}_${pvalue}.txt \
-o ${outDir}/pval_${pvalue}

mv ${outDir}/pval_${pvalue}/${name}.bed ${outDir}/pval_${pvalue}/${name}_${pvalue}.bed 

# Create a bedgraph from the bed output
echo $(date +"[%b %d %H:%M:%S] Create a bedgraph from bed output...")

bedtools sort \
-i ${outDir}/pval_${pvalue}/${name}_${pvalue}.bed \
| bedtools merge \
-i - \
-o max \
-c 5 \
> ${outDir}/pval_${pvalue}/${name}_${pvalue}.bdg

# Create a bigwig from the bedgraph output
echo $(date +"[%b %d %H:%M:%S] Create a bigwig from bedgraph output...")

/Users/CL_Shared/programs/ucsc/bedGraphToBigWig \
${outDir}/pval_${pvalue}/${name}_${pvalue}.bdg \
${chromSizes} \
${outDir}/pval_${pvalue}/${name}_${pvalue}.bw

echo $(date +"[%b %d %H:%M:%S] Done")
