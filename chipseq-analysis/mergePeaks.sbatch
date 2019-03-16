#!/bin/bash

## Script for merging peaks from all samples in a batch
## Date: 16 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/macs2/filtered_peaks \
## projectName=project20xx_chipseq\
## sbatch mergePeaks.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J mergePeaks
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load bedtools

# Merge peaks from all samples
pwd; hostname; date

echo "Merging peaks from all samples..."

# Go to dir
cd ${inDir}

# Concatenate all narrowPeak regions from all samples
cat *.narrowPeak.filtered > "$projectName"_all.tmp

# Sort and merge intervals within 100bp from each other
bedtools sort -i "$projectName"_all.tmp > "$projectName"_all_sorted.tmp
bedtools merge -i "$projectName"_all_sorted.tmp -d 100 > "$projectName"_all_sorted_merged100bp.tmp

# Label the regions
cat "$projectName"_all_sorted_merged100bp.tmp \
| awk '{print $0 "\t" "region"NR}' \
> all_"$projectName".narrowPeaks

# Remove tmp files
rm *.tmp

echo $(date +"[%b %d %H:%M:%S] Done")