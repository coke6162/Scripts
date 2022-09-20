#!/bin/bash

## Script for running fimo
## Date: September 2 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/genome-setup
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/fa/bosTau9.fa \
## chromSizes=/Shares/CL_Shared/db/genomes/bosTau9/fa/bosTau9.chrom.sizes \
## pvalue=1e-4 \
## motifFile=/Users/coke6162/repos/scripts/genome-setup/MA0137.3.meme \
## name=MA0137.3 \
## outDir=/Shares/CL_Shared/db/genomes/bosTau9/fimo/MA0137.3 \
## sbatch fimo.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

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

# Make output directories
pwd; hostname; date

echo "fimo version: "$(fimo --version)
echo "bedtools version: "$(bedtools --version)

echo "Making output directories..."

mkdir -p ${outDir}/${name}/pval_${pvalue}

# Run fimo
echo $(date +"[%b %d %H:%M:%S] Running fimo...")

fimo \
--text \
--parse-genomic-coord \
--thresh ${pvalue} \
${motifFile} \
${fastaFile} \
> ${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.txt

# Convert output to bed & add pval as a file suffix
echo $(date +"[%b %d %H:%M:%S] Convert fimo output to bed format...")

mkdir -p ${outDir}/${name}/pval_${pvalue}

/opt/python/2.7.14/bin/python2.7 /Shares/CL_Shared/scripts/ed/fimo_analysis.py \
-i ${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.txt \
-o ${outDir}/${name}/pval_${pvalue}

mv ${outDir}/${name}/pval_${pvalue}/${name}.bed ${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.bed 

# Create a bedgraph from the bed output
echo $(date +"[%b %d %H:%M:%S] Create a bedgraph from bed output...")

bedtools sort \
-i ${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.bed \
| bedtools merge \
-i - \
-o max \
-c 5 \
> ${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.bdg

# Create a bigwig from the bedgraph output
echo $(date +"[%b %d %H:%M:%S] Create a bigwig from bedgraph output...")

/Shares/CL_Shared/programs/ucsc/bedGraphToBigWig \
${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.bdg \
${chromSizes} \
${outDir}/${name}/pval_${pvalue}/${name}_${pvalue}.bw

echo $(date +"[%b %d %H:%M:%S] Done")
