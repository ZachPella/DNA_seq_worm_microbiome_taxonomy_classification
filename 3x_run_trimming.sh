#!/bin/bash
#SBATCH --job-name=trim_unmapped
#SBATCH --output=trim_%A_%a.out
#SBATCH --error=trim_%A_%a.err
#SBATCH --array=1-9
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4

cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

module load trimmomatic
module load java


# Get list of R1 files
R1_FILES=($(ls fastq_for_assembly/*.unmapped.R1.fastq))

# Get current file based on array task ID
INDEX=$((SLURM_ARRAY_TASK_ID-1))
R1_FILE=${R1_FILES[$INDEX]}
R2_FILE=${R1_FILE/R1.fastq/R2.fastq}

# Extract sample name
SAMPLE=$(basename ${R1_FILE} .unmapped.R1.fastq)

# Create output directory
mkdir -p trimmed_fastq

# Run Trimmomatic (without adapter trimming)
java -jar $TM_HOME/trimmomatic.jar PE -threads 4 \
    ${R1_FILE} ${R2_FILE} \
    trimmed_fastq/${SAMPLE}.trimmed.R1.fastq trimmed_fastq/${SAMPLE}.unpaired.R1.fastq \
    trimmed_fastq/${SAMPLE}.trimmed.R2.fastq trimmed_fastq/${SAMPLE}.unpaired.R2.fastq \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
