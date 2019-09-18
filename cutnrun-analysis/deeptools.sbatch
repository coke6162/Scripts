#!/bin/bash

## Script for running deeptools 
## Date: 18 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## fileName=MER41_BT \
## bedFile=/Users/CL_Shared/db/genomes/bosTau9/repeats/eachRep/MER41_BT.bed \
## wigFile=bigwigs_all.txt \
## refDir=/Users/coke6162/repos/project_files/20190820_CUTnRUN_MDBK \
## outDir=/Users/coke6162/project20xx_cutnrun/deeptools \
## sbatch deeptools.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J deeptools
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# load modules
module load singularity

# define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
wigs=$(cat ${refDir}/${wigFile})
samplesLabel=$(cat ${refDir}/${wigFile} | xargs -n 1 basename | sed 's/.bw//g')

# run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate data underlying heatmap
windowLeft=4000
windowRight=4000
binSize=10
numCPU=16

singularity exec --bind /Users/CL_Shared ${deeptools} \
computeMatrix reference-point \
--referencePoint TSS \
--scoreFileName ${wigs} \
--regionsFileName ${bedFile} \
--beforeRegionStartLength ${windowLeft} \
--afterRegionStartLength ${windowRight} \
--binSize ${binSize} \
--missingDataAsZero \
-o ${outDir}/${fileName}.mat.gz \
-p ${numCPU}

echo $(date +"[%b %d %H:%M:%S] Plot heatmap...")

# Use "plotHeatmap" to create a png or pdf
zMin=0
yMin=0

singularity exec --bind /Users/CL_Shared ${deeptools} \
plotHeatmap \
-m ${outDir}/${fileName}.mat.gz \
--outFileName ${outDir}/${fileName}.png \
--outFileSortedRegions ${outDir}/${fileName}.dt.bed \
--outFileNameMatrix ${outDir}/${fileName}.matrix.tab \
--sortRegions descend \
--colorMap Blues \
--zMin ${zMin} \
--yMin ${yMin} \
--regionsLabel ${fileName}.bed \
--samplesLabel ${samplesLabel} \
--plotTitle ${fileName}_heatmap \

echo $(date +"[%b %d %H:%M:%S] Done")
