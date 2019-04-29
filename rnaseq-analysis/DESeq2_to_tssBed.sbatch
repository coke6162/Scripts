#!/bin/bash

## Script for getting tss coordinates from DESeq2 output
## Date: 29 April 2019

## Example usage:
## inDir=/Users/coke6162/project20xx_rnaseq/DESeq2 \
## outDir=/Users/coke6162/project20xx_rnaseq/DESeq2/bed \
## bedFile=/Users/CL_Shared/db/genomes/mm10/annotations/gencode.vM17.annotation.tss.bed \
## sbatch --array 0-3 DESeq2_to_tssBed.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=0:30:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J fastqc
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*_padj.txt | xargs -n 1 basename | sed --expression='s/_padj.txt//g'))

# Get tss coordinates for differentially expressed genes
pwd; hostname; date

echo "Processing files: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Generating tss bed files from DESeq2 output...")

awk '(NR>1){print $1}' ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_padj.txt \
| awk 'FNR==NR {a[$4]=$0; next}; $1 in a {print a[$0]}' \
${bedFile} - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_padj_tss.bed

awk '(NR>1){print $1}' ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_padj_up.txt \
| awk 'FNR==NR {a[$4]=$0; next}; $1 in a {print a[$0]}' \
${bedFile} - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_padj_up_tss.bed

awk '(NR>1){print $1}' ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_padj_down.txt \
| awk 'FNR==NR {a[$4]=$0; next}; $1 in a {print a[$0]}' \
${bedFile} - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_padj_down_tss.bed

echo $(date +"[%b %d %H:%M:%S] Done")