#!/bin/bash

## Script to search for repeat enrichment using giggle
## Date: 16 September 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_cutnrun/macs2 \
## outDir=/Users/coke6162/project20xx_cutnrun/giggle/non-DE/allRep \
## repeatIndexDir=/Shares/CL_Shared/db/giggle/bosTau8/repeats/indexed \
## repeatPath=/Users/CL_Shared/db/giggle/bosTau8/repeats/sorted/ \
## genomeSize=2660906405 \
## suffix=narrowPeak \
## sbatch --array 0-3 giggleSearch.sbatch

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
queries=($(ls ${inDir}/*${suffix} | xargs -n 1 basename))

# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "samtools version: "$(samtools --version)
echo "Processing peak file: "${queries[$SLURM_ARRAY_TASK_ID]}

# make a temporary directory
echo $(date +"[%b %d %H:%M:%S] Making a temporary directory...")

mkdir -p ${inDir}/tmp

# add '.bed' suffix
echo $(date +"[%b %d %H:%M:%S] Adding '.bed' suffix...")

cp ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} ${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]}.bed

# bgzip each bed file
echo $(date +"[%b %d %H:%M:%S] Starting bgzip...")

bgzip ${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]}.bed

# run giggle search
echo $(date +"[%b %d %H:%M:%S] Starting giggle search...")

/Shares/CL_Shared/scripts/atma/giggle/bin/giggle search \
-i ${repeatIndexDir} \
-q ${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]}.bed.gz \
-s -g ${genomeSize} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.tmp

# remove full index name from output
echo $(date +"[%b %d %H:%M:%S] Removing full repeat name...")

awk '(NR>1){print $0}' ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.tmp \
| sed 's/.bed.gz//g' \
| sed "s|$repeatPath||g" \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats

# sort output by descending giggle score
# note: only use unsorted output as input into makeGiggleMatrix.sbatch!
echo $(date +"[%b %d %H:%M:%S] Sorting by descending combo score...")

sort -k8,8rn ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.sorted

# remove temporary file
echo $(date +"[%b %d %H:%M:%S] Removing intermediate file...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.tmp

# remove '.bed.gz' files
echo $(date +"[%b %d %H:%M:%S] Removing '.bed.gz' files...")

rm ${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]}.bed.gz

echo $(date +"[%b %d %H:%M:%S] Done")