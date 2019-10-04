#!/bin/bash

## Script for running deeptools to compute tss enrichment matrices
## Date: 3 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## fileName=bosTau8_ensGene_tss \
## wigDir=/Users/coke6162/project20xx_atacseq/bigwigs \
## bedFile=/Users/CL_Shared/db/genomes/bosTau8/annotation/bosTau8.ensGene.tss.bed \
## outDir=/Users/coke6162/project20xx_atacseq/bowtie2/tss_enrich \
## sbatch --array 0-3 deeptools_tssEnrich.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J deeptools
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# load modules
module load singularity

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# Define query files
queries=($(ls ${wigDir}/*_treat_pileup.bw | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate enrichment matrix
windowLeft=2000
windowRight=2000
binSize=10
numCPU=16

singularity exec --bind /Users/CL_Shared ${deeptools} \
computeMatrix reference-point \
--referencePoint TSS \
--scoreFileName ${wigDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--regionsFileName ${bedFile} \
--outFileName ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_treat_pileup.bw}_${fileName}.matrix.gz \
--outFileNameMatrix ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_treat_pileup.bw}_${fileName}.matrix \
--beforeRegionStartLength ${windowLeft} \
--afterRegionStartLength ${windowRight} \
--binSize ${binSize} \
--missingDataAsZero \
--numberOfProcessors ${numCPU} \

echo $(date +"[%b %d %H:%M:%S] Plot heatmap...")