#!/bin/bash
#SBATCH -p short                   # Partition or queue
#SBATCH --job-name=project20xx_rnaseq_fastqc_trim         # Job name
#SBATCH --mail-type=END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=coke6162@colorado.edu     # Where to send mail
#SBATCH --ntasks=1                 # Run on a single CPU
#SBATCH --cpus-per-task=8         # cpus
#SBATCH --mem=4gb                # Memory limit
#SBATCH --time=1:00:00            # Time limit hrs:min:sec
#SBATCH --output=/Users/coke6162/sbatch_script/output/project20xx_rnaseq_fastqc_trim_%A_%a.out      # Standard output log
#SBATCH --error=/Users/coke6162/sbatch_script/error/project20xx_rnaseq_fastqc_trim_%A_%a.err       # Standard error log

# Load FastQC
module load fastqc

# Define constant variables
numThreads=8

# Define query files
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename))

# Run FastQC
pwd; hostname; date

echo "Running fastqc..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

fastqc -o ${outDir} -f fastq -t ${numThreads} ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

date
echo "Done"

# Command line entry
### INDIR=/Users/coke6162/project20xx/trimmomatic \
### OUTDIR=/Users/coke6162/project20xx/trimmomatic/fastqc \
### sbatch --array 0-103 /Users/coke6162/project20xx/scripts/fastqc_trim.sbatch
