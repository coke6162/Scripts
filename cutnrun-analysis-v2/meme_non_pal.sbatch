#!/bin/bash

## Script for running meme in non-palindromic mode (default)
## Date: 15 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/input/fa \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/non_pal \
## dbFile=/Shares/CL_Shared/programs/meme-5.1.0/db/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme \
## sbatch --array=0-0 meme_non_pal.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=6:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J meme
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.fa | xargs -n 1 basename | sed 's/.fa//g' | uniq))

# make output directory
pwd; hostname; date

echo "Making output & summary directories..."

mkdir -p ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/tomtom
mkdir -p ${outDir}/summary

# run meme
echo $(date +"[%b %d %H:%M:%S] Running meme...")

echo "meme version: "$(/Shares/CL_Shared/bin/meme/bin/meme -version)
echo "tomtom version: "$(/Shares/CL_Shared/bin/meme/bin/tomtom -version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa

/Shares/CL_Shared/bin/meme/bin/meme \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa \
-oc ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-dna \
-mod zoops \
-minw 6 \
-maxw 20 \
-nmotifs 20 \
-markov_order 1

# run tomtom
echo $(date +"[%b %d %H:%M:%S] Running tomtom...")

echo "Searching against database: "${dbFile}

/Shares/CL_Shared/bin/meme/bin/tomtom \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/meme.txt \
${dbFile} \
-oc ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/tomtom

# copy results to convenient summary dir
echo $(date +"[%b %d %H:%M:%S] Copying results to convenient summary directory...")

cp ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/meme.html ${outDir}/summary/${queries[$SLURM_ARRAY_TASK_ID]}_meme.html
cp ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/tomtom/tomtom.html ${outDir}/summary/${queries[$SLURM_ARRAY_TASK_ID]}_tomtom.html

echo $(date +"[%b %d %H:%M:%S] Done")