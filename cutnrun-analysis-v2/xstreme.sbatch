#!/bin/bash

## Script for running xstreme
## Date: 18 September 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/input/fa \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/11_meme_output/non_pal \
## dbFile=/Users/coke6162/programs/meme-5.4.1/db/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme \
## sbatch --array=0-0 xstreme.sbatch

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=4-0:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J xstreme
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

echo "Making summary directory..."

mkdir -p ${outDir}/summary

# run xstreme
echo $(date +"[%b %d %H:%M:%S] Running streme...")

echo "xstreme version: "$(streme -version)
echo "tomtom version: "$(tomtom -version)

echo "Processing file: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa

xstreme \
--verbosity 2 \
--oc ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--dna \
--p ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.fa \
--m ${dbFile} \
--minw 6 \
--maxw 20 \
--streme-nmotifs 20 \
--align center

# copy results to convenient summary dir
echo $(date +"[%b %d %H:%M:%S] Copying results to convenient summary directory...")

cp ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}/xstreme.html ${outDir}/summary/${queries[$SLURM_ARRAY_TASK_ID]}_xstreme.html

echo $(date +"[%b %d %H:%M:%S] Done")
