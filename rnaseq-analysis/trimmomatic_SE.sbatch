#!/bin/bash
#SBATCH -p short                   # Partition or queue
#SBATCH --job-name=project20xx_rnaseq_trim         # Job name
#SBATCH --mail-type=END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=coke6162@colorado.edu     # Where to send mail
#SBATCH --ntasks=1                 # Run on a single CPU
#SBATCH --cpus-per-task=8         # cpus
#SBATCH --mem=32gb                # Memory limit
#SBATCH --time=1:00:00            # Time limit hrs:min:sec
#SBATCH --output=/Users/coke6162/sbatch_script/output/project20xx_rnaseq_trim_%A_%a.out      # Standard output log
#SBATCH --error=/Users/coke6162/sbatch_script/error/project20xx_rnaseq_trim_%A_%a.err       # Standard error log

# Set constant variables
trimExe=/opt/trimmomatic/0.36/trimmomatic-0.36.jar
numThreads=8
adapterFile=TruSeq3-PE.fa # NexteraPE-PE.fa, TruSeq2-PE.fa, TruSeq3-PE-2.fa, TruSeq3-PE.fa

# Set query files
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename))

# Run Trimmomatic
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Starting trimmomatic..."

java -jar ${trimExe} PE -threads ${numThreads} \
-trimlog ${logOutDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}_trimlog.txt \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}_trimmed.fastq.gz \
ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/${adapterFile}:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

echo $(date +"[%b %d %H:%M:%S] Done")

# Command line entry
### INDIR=/Users/coke6162/project20xx_rnaseq/fastq \
### OUTDIR=/Users/coke6162/project20xx_rnaseq/trimmomatic \
### LOGOUTDIR=/Users/coke6162/project20xx_rnaseq/trimmomatic/log \
### sbatch --array 0-8 /Users/coke6162/project20xx_rnaseq/scripts/trim.sbatch
