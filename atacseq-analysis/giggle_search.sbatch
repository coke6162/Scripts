#!/bin/bash

## Script to search for repeat enrichment using giggle
## Date: 8 May 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_atacseq/macs2 \
## outDir=/Users/coke6162/project20xx_atacseq/giggle/allRep \
## repeatIndexDir=/Users/CL_Shared/db/giggle/bosTau8/repeats/indexed \
## repeatPath=/Users/CL_Shared/db/giggle/bosTau8/repeats/sorted/ \
## sbatch --array 0-3 giggle_search.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --propagate=NONE

# Job name and output
#SBATCH -J giggle_search
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# define query files
queries=($(ls ${inDir}/*.narrowPeak | xargs -n 1 basename))

# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "samtools version: "$(samtools --version)
echo "Processing peak file: "${queries[$SLURM_ARRAY_TASK_ID]}

# extract first three columns from peak file to make a simple bed
echo $(date +"[%b %d %H:%M:%S] Converting peak file to simple bed...")

cat ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| awk '{print $1 "\t" $2 "\t" $3}' \
> ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed

# bgzip each bed file
echo $(date +"[%b %d %H:%M:%S] Starting bgzip...")

bgzip ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed

# run giggle search
echo $(date +"[%b %d %H:%M:%S] Starting giggle search...")

giggle search \
-i ${repeatIndexDir} \
-q ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed.gz \
-s \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bed.gz}.giggleStats.tmp

# remove full index name from output & sort by combo score
echo $(date +"[%b %d %H:%M:%S] Removing full repeat name and sort by combo score...")

cat ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bed.gz}.giggleStats.tmp \
| sed 's/.bed.gz//g' \
| sed "s|$repeatPath||g" \
| sort -k8,8rn \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_slop.bed.gz}.giggleStats

# remove temporary file
echo $(date +"[%b %d %H:%M:%S] Removing intermediate file...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bed.gz}.giggleStats.tmp

echo $(date +"[%b %d %H:%M:%S] Done")