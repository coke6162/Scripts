#!/bin/bash

## Script for running bowtie2 index
## Date: 1 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## indexDir=/Users/CL_Shared/db/genomes/bosTau8/index \
## inDir=/Users/coke6162/project20xx_atacseq/bbduk \
## outDir=/Users/coke6162/project20xx_atacseq/bowtie2 \
## name=bosTau8 \
## sbatch --array=0-3 bowtie2.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J bowtie2
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Module load bowtie2 and samtools
module load bowtie
module load samtools

# Define query files
queries=($(ls ${inDir}/*fastq.gz | xargs -n 1 basename | sed --expression='s/_R1_paired.fastq.gz//g' | sed --expression='s/_R2_paired.fastq.gz//g' | uniq))

# Align reads to bosTau8
pwd; hostname; date

echo "Aligning reads to the genome & outputting sam files"
echo "Bowtie2 version: "$(bowtie2 --version)
echo "Samtools version: "$(samtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bowtie2 --end-to-end --very-sensitive -X 1000 --fr -x ${indexDir}/${name} \
-1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz \
-2 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

# Convert sam files to bam files, filtering out low quality (<10) & unmapped (through '-F 4' option) reads
echo "Convert sams to bams"
samtools view -@ 8 -Sb -q 10 -F 4 ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bam

# Sort bam files
echo "Sort bams"
samtools sort -@ 8 ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Generate flagstat files from sorted bam files
echo "Make flagstat files"
samtools flagstat ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_flagstat.txt

# Index sorted bam files
echo "Index sorted bams"
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Determine percentage of total reads that mapped to mitochondrial genes
echo "Calculate percentage of total reads mapping to mitochondrial genes"
chrMreads=`samtools view -c ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam chrM`
totalReads=`samtools view -c ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam`
fractionMreads=`echo "100 * ${chrMreads} / ${totalReads}" | bc -l`
touch ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${queries[$SLURM_ARRAY_TASK_ID]} >> ${queries[$SLURM_ARRAY_TASK_ID]}/_chrMreadsFraction.txt
echo ${totalReads} 'total mapped reads' >> ${queries[$SLURM_ARRAY_TASK_ID]}/_chrMreadsFraction.txt
echo ${chrMreads} 'mitochondrial reads' >> ${queries[$SLURM_ARRAY_TASK_ID]}/_chrMreadsFraction.txt
echo ${fractionMreads} 'percentage of mitochondrial reads from total mapped reads' >> ${queries[$SLURM_ARRAY_TASK_ID]}/_chrMreadsFraction.txt

echo $(date +"[%b %d %H:%M:%S] Done")

### Explanation of bowtie2 options
### '--end-to-end' - specifies end-to-end alignment (options: '--end-to-end' and '--local')
### '--sensitive' - is a preset that specifies subset of settings for "sensitive" end-to-end alignment (options: '--very-fast', '--fast', '--sensitive', '--very-sensitive')
### '-X <int>' - specifies maximum fragment length of valid paired-end alignments; allows for gaps between pairs
### '--fr' - specifies upstream/downstream mate orientations for a valid paired-end alignment against the forward reference strand; e.g. if mate 1 appears upstream of reverse complement of mate 2, the alignment is valid
### '-x <bt2-index>' - used to specify the basename of the index for the reference genome; the basename is the name of any index files up to but not including the final '.1.bt2'
### '-1 <r1>' and '-2 <r2>' - used to specify paired input

### Explanation of samtools options
### '-@ <int>' - number of threads
### '-Sb' - convert sam to bam
### '-q <int>' - quality threshold (filter out reads with quality below int)
### '-F 4' - exclude unmapped reads
### '-c' - print only the count of matching records