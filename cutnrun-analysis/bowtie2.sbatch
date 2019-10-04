#!/bin/bash

## Script for running bowtie2 
## Date: 11 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## indexDir=/Users/CL_Shared/db/genomes/bosTau9/index \
## inDir=/Users/coke6162/project20xx_cutnrun/bbduk \
## outDir=/Users/coke6162/project20xx_cutnrun/bowtie2 \
## name=bosTau9.main \
## sbatch --array=0-3 bowtie2.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
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
module load bedtools

# Define query files
queries=($(ls ${inDir}/*fastq.gz | xargs -n 1 basename | sed --expression='s/_R1_paired.fastq.gz//g' | sed --expression='s/_R2_paired.fastq.gz//g' | uniq))

# Align reads to bosTau8
pwd; hostname; date

echo "Aligning reads to the genome & outputting sam files"
echo "Bowtie2 version: "$(bowtie2 --version)
echo "Samtools version: "$(samtools --version)
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bowtie2 \
--local \
--very-sensitive-local \
--no-unal \
--no-mixed \
--no-discordant \
-I 10 \
-X 700 \
--fr \
--threads ${numThreads} \
-x ${indexDir}/${name} \
-1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz \
-2 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sam

# Convert sam files to bam files, filtering out low quality (<10) & unmapped (through '-F 4' option) reads
echo "Convert sams to bams"
samtools view \
-@ ${numThreads} \
-Sb \
-q 10 \
-F 4 \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.bam

# Sort bam files
echo "Sort bams"
samtools sort \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam

# Generate flagstat files from sorted bam files
echo "Make flagstat files"
samtools flagstat \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_flagstat.txt

# Index sorted bam files
echo "Index sorted bams"
samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam.bai

# Convert to BEDPE & filter out fragments greater than 1kb
echo "Converting bam to bed"

samtools sort \
-n ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam \
| samtools fixmate - - \
| bamToBed -bedpe -i - \
| awk '$1==$4 && $6-$2 < 1000 {print $0}' - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_d1000.sorted.bed

# Select the 5' and 3' coordinates of the read pair to generate a new BED3 file
cut -f 1,2,6 ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_d1000.sorted.bed \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_d1000_fragments.sorted.bed

# Determine percentage of total reads that mapped to mitochondrial genes
echo "Calculate percentage of total reads mapping to mitochondrial genes"
chrMreads=`samtools view -c ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam chrM`
totalReads=`samtools view -c ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}.sorted.bam`
fractionMreads=`echo "100 * ${chrMreads} / ${totalReads}" | bc -l`
touch ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_chrMreadsFraction.txt
echo ${queries[$SLURM_ARRAY_TASK_ID]} >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_chrMreadsFraction.txt
echo ${totalReads} 'total mapped reads' >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_chrMreadsFraction.txt
echo ${chrMreads} 'mitochondrial reads' >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_chrMreadsFraction.txt
echo ${fractionMreads} 'percentage of mitochondrial reads from total mapped reads' >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_${name}_chrMreadsFraction.txt

echo $(date +"[%b %d %H:%M:%S] Done")

### Explanation of bowtie2 options
### '--lcaol' - specifies local alignment (options: '--end-to-end' and '--local')
### '--very-sensitive-local' - is a preset that specifies subset of settings for "very sensitive" local alignment (options: '--very-fast-local', '--fast-local', '--sensitive-local', '--very-sensitive-local')
### '--no-unal' - suppress SAM records for reads that fail to align
### '--no-mixed' - if bowtie2 cannot find a concordant or discordant alignment for a pair, do not try to align individual mates
### '--no-discordant' - disables discordant mapping (i.e. when mates map uniquely but in a way that does not satisfy paired-end constraints, i.e. --fr/--rf/--ff, -I, -X)
### '--I' - specified minimum fragment length of valid paired-end alignments
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