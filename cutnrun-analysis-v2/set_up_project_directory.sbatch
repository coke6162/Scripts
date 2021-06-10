#!/bin/bash

## Script to set up directory layout on /scratch
## Date: 9 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## projectDir=/scratch/Users/coke6162/yyyymmdd_projectName \
## sbatch set_up_project_directory.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=0:5:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J setup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# run the thing
pwd; hostname; date

echo "Project: "${projectDir}
echo $(date +"[%b %d %H:%M:%S] Making project directories...")

# make a new project dir
mkdir -p ${projectDir}

# make necessary subdirs
mkdir -p ${projectDir}/1_raw_fastq
mkdir -p ${projectDir}/1_raw_fastq/control
mkdir -p ${projectDir}/1_raw_fastq/histone
mkdir -p ${projectDir}/1_raw_fastq/tf

mkdir -p ${projectDir}/2_trimmed_fastq
mkdir -p ${projectDir}/2_trimmed_fastq/control
mkdir -p ${projectDir}/2_trimmed_fastq/histone
mkdir -p ${projectDir}/2_trimmed_fastq/tf

mkdir -p ${projectDir}/3_bams
mkdir -p ${projectDir}/3_bams/control
mkdir -p ${projectDir}/3_bams/histone
mkdir -p ${projectDir}/3_bams/tf

mkdir -p ${projectDir}/4_macs2_output_PE
mkdir -p ${projectDir}/4_macs2_output_PE/histone
mkdir -p ${projectDir}/4_macs2_output_PE/tf

mkdir -p ${projectDir}/5_macs2_output_SE
mkdir -p ${projectDir}/5_macs2_output_SE/histone
mkdir -p ${projectDir}/5_macs2_output_SE/tf

mkdir -p ${projectDir}/6_macs2_output_merged
mkdir -p ${projectDir}/6_macs2_output_merged/histone
mkdir -p ${projectDir}/6_macs2_output_merged/tf

mkdir -p ${projectDir}/7_bigwigs

mkdir -p ${projectDir}/8_seacr_output
mkdir -p ${projectDir}/8_seacr_output/control
mkdir -p ${projectDir}/8_seacr_output/histone
mkdir -p ${projectDir}/8_seacr_output/tf

mkdir -p ${projectDir}/9_signal_heatmaps

mkdir -p ${projectDir}/10_giggle_output

mkdir -p ${projectDir}/reports
mkdir -p ${projectDir}/reports/raw_fastqc
mkdir -p ${projectDir}/reports/raw_multiqc
mkdir -p ${projectDir}/reports/trimmed_fastqc
mkdir -p ${projectDir}/reports/trimmed_multiqc
mkdir -p ${projectDir}/reports/bam_fragment_size
mkdir -p ${projectDir}/reports/frip_scores

echo $(date +"[%b %d %H:%M:%S] Done")