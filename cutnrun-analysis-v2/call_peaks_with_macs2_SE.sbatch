#!/bin/bash

## Script for calling peaks with macs2 in SE mode (against IgG control)
## Date: 11 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/3_bams/histone \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/4_macs2_output_SE/histone \
## controlFile=/scratch/Users/coke6162/yyyymmdd_projectName/3_bams/control/*.bam \
## genomeSize=2.72e9 \
## sbatch --array=0-0 call_peaks_with_macs2_SE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J macs2_SE
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load python/2.7.14/MACS/2.1.1

# Define query files
queries=($(ls ${inDir}/*_filtered.sorted.bam | xargs -n 1 basename | sed 's/_filtered.sorted.bam//g' | sort | uniq))

# Run the thing
pwd; hostname; date

echo "macs2 version: "$(macs2 --version)

echo "Target file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Control file: "${controlFile}

echo "Starting macs2..."

macs2 callpeak \
--format BAM \
--shift=-75 \
--extsize=150 \
--treatment ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--control ${controlFile} \
--pvalue 0.01 \
--name ${queries[$SLURM_ARRAY_TASK_ID]} \
--outdir ${outDir} \
--gsize ${genomeSize} \
--bdg \
--SPMR \
--call-summits

# Remove files
echo "[%b %d %H:%M:%S] Removing files..."

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_peaks.xls
#rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_model.r
#rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_control_lambda.bdg

echo $(date +"[%b %d %H:%M:%S] Done")

### Explanation of macs2 options
### '--format BAM' - specifies format of input file
### '--shift=<int>' - shift read ends by INT; if negative, shift 3' -> 5'
### '--extsize=<int>' - extend read size by INT
### '--treatment <file>' - pulldown file
### '--control <file>' - control file; for CUT&RUN this is IgG
### '--pvalue <int>' - set pvalue threshold of <int>; default is qval 0.05
### '--name <string>' - prefix for output files
### '--outdir <dir>' - path to output directory
### '--gsize <int>' - approximate genome size; default is 'hs' or 2.7e9
### '--bdg' - write bedgraph file
### '--SPMR' - normalize bedgraph by signal per million reads
### '--call-sumits' - find subpeak summits for each peak