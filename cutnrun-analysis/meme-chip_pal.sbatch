#!/bin/bash

## Script for running meme-chip
## Date: 20 December 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/yyyymmdd_CUTnRUN_projectName/seacr/stringent/slopSummits/b250/fa \
## outDir=/Users/coke6162/yyyymmdd_CUTnRUN_projectName/meme/meme-chip/seacr/stringent/slopSummits_b250 \
## dbFile=/Shares/CL_Shared/programs/meme-5.1.0/db/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme \
## sbatch --array=0-5 meme-chip_pal.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=5-0:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J meme-chip
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.fa | xargs -n 1 basename | sed 's/.fa//g' | uniq))

# make output directory
pwd; hostname; date

echo "Making output & summary directories..."

mkdir -p ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}
mkdir -p ${outDir}/summary

# run the thing
echo $(date +"[%b %d %H:%M:%S] Starting meme-chip...")
echo "meme-chip version: "$(meme-chip -version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa
echo "Searching against database: "${dbFile}

/Shares/CL_Shared/bin/meme/bin/meme-chip \
-oc ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-ccut 200 \
-order 1 \
-db ${dbFile} \
-meme-mod zoops \
-meme-minw 8 \
-meme-maxw 14 \
-meme-nmotifs 10 \
-meme-pal \
-spamo-skip \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa

# copy results to convenient summary dir
echo $(date +"[%b %d %H:%M:%S] Copying results to convenient summary directory...")

cp ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/meme-chip.html ${outDir}/summary/${queries[$SLURM_ARRAY_TASK_ID]}_meme-chip.html

echo $(date +"[%b %d %H:%M:%S] Done")