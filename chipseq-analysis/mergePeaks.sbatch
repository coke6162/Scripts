#!/bin/bash

## Script for merging peaks from all samples in a batch
## Date: 11 July 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/chipseq-analysis
## pulldown=STAT1 \
## treatment=1.5h_IFNg \
## fileType=narrowPeak.filtered \
## inDir=/Users/coke6162/project20xx_chipseq/macs2 \
## outDir=/Users/coke6162/project20xx_chipseq/DESeq2 \
## sbatch mergePeaks.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J mergePeaks
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load bedtools

# Concatenate all narrowPeak regions from pol2 pulldowns
pwd; hostname; date

echo "Merging narrowPeak regions from all samples"
echo "bedtools version: "$(bedtools --version)

cat ${inDir}/*${treatment}*${pulldown}*${fileType} > ${outDir}/all_${pulldown}_${treatment}_${fileType}.tmp

echo $(date +"[%b %d %H:%M:%S]")

# Sort
echo "Sorting regions..."

bedtools sort \
-i ${outDir}/all_${pulldown}_${treatment}_${fileType}.tmp \
> ${outDir}/all_${pulldown}_${treatment}_${fileType}_sorted.tmp

echo $(date +"[%b %d %H:%M:%S]")

# Merge intervals within 100bp from each other
echo "Merging intervals within 100bp from each other..."

bedtools merge \
-i ${outDir}/all_${pulldown}_${treatment}_${fileType}_sorted.tmp \
-d 100 \
> ${outDir}/all_${pulldown}_${treatment}_${fileType}_sorted_merged100bp.tmp

echo $(date +"[%b %d %H:%M:%S]")

# Label the regions
echo "Labeling the regions..."

cat ${outDir}/all_${pulldown}_${treatment}_${fileType}_sorted_merged100bp.tmp \
| awk '{print $0 "\t" "region"NR}' \
> ${outDir}/all_${pulldown}_${treatment}_${fileType}

echo $(date +"[%b %d %H:%M:%S]")

# Remove tmp files
echo "Removing intermediate files..."

rm ${outDir}/all_${pulldown}_${treatment}_${fileType}.tmp
rm ${outDir}/all_${pulldown}_${treatment}_${fileType}_sorted.tmp
rm ${outDir}/all_${pulldown}_${treatment}_${fileType}_sorted_merged100bp.tmp

echo $(date +"[%b %d %H:%M:%S] Done")