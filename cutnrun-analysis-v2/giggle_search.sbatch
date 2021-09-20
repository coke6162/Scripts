#!/bin/bash

## Script to search for repeat enrichment using giggle
## Date: 14 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/6_macs2_output_merged/histone \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/10_giggle_output \
## repeatIndexDir=/Shares/CL_Shared/db/giggle/bosTau9/repeats/indexed \
## repeatPath=/Shares/CL_Shared/db/giggle/bosTau9/repeats/sorted/ \
## genomeSize=2715853792 \
## suffix="_mergedpeaks.bed" \
## sbatch --array 0-0 giggle_search.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=5MB
#SBATCH --propagate=NONE

# Job name and output
#SBATCH -J giggle_search
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# important note: make sure the giggle executable is in your PATH

# load modules
module load samtools

# define query files
queries=($(ls ${inDir}/*${suffix} | xargs -n 1 basename))

# make a temporary directory
pwd; hostname; date

echo "Making a temporary directory..."

mkdir -p ${inDir}/tmp
cd ${inDir}/tmp

# add '.bed' suffix
echo $(date +"[%b %d %H:%M:%S] Adding '.bed' suffix...")

cp ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} ${inDir}/tmp/${queries[$SLURM_ARRAY_TASK_ID]}.bed

# bgzip each bed file
echo $(date +"[%b %d %H:%M:%S] Compressing files ending in '.bed'...")

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

echo $(date +"[%b %d %H:%M:%S] Done!")