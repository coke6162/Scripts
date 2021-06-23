#!/bin/bash

## Script for running streme in non-palindromic mode (default)
## Date: 22 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/input/fa \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/non_pal \
## dbFile=/Shares/CL_Shared/programs/streme-5.1.0/db/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.streme \
## sbatch --array=0-0 streme.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=6:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J streme
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

# run streme
echo $(date +"[%b %d %H:%M:%S] Running streme...")

echo "streme version: "$(streme -version)
echo "tomtom version: "$(tomtom -version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa

streme \
--verbosity 2 \
--oc ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--dna \
--p ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa \
--minw 6 \
--maxw 20 \
--nmotifs 20 \
--totallength 4000000

# run tomtom
echo $(date +"[%b %d %H:%M:%S] Running tomtom...")

echo "Searching against database: "${dbFile}

tomtom \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/streme.txt \
${dbFile} \
-oc ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/tomtom

# copy results to convenient summary dir
echo $(date +"[%b %d %H:%M:%S] Copying results to convenient summary directory...")

cp ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/streme.html ${outDir}/summary/${queries[$SLURM_ARRAY_TASK_ID]}_streme.html
cp ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/tomtom/tomtom.html ${outDir}/summary/${queries[$SLURM_ARRAY_TASK_ID]}_streme_tomtom.html

echo $(date +"[%b %d %H:%M:%S] Done")