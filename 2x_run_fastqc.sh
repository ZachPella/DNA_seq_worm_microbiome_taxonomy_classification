#!/bin/bash
#SBATCH --job-name=fastqc_unmapped
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --error=fastqc_%A_%a.err
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-18

cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

module load fastqc

# Create output directory
mkdir -p fastqc_results

# Get list of non-empty FASTQ files (exclude singleton files which are empty)
FASTQ_FILES=($(find fastq_for_assembly/ -name "*.fastq" -size +0c | sort))

# Calculate which file this array task should process
FILE_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
CURRENT_FILE=${FASTQ_FILES[$FILE_INDEX]}

# Run FastQC on the specific file for this array task
if [ -f "$CURRENT_FILE" ]; then
    echo "Processing file: $CURRENT_FILE"
    fastqc -t 4 -o fastqc_results "$CURRENT_FILE"
else
    echo "File not found: $CURRENT_FILE"
    exit 1
fi

# Only run MultiQC once, when the last array job completes
if [ $SLURM_ARRAY_TASK_ID -eq 18 ]; then
    # Wait a bit to ensure other jobs have finished
    sleep 30

    module load multiqc
    multiqc fastqc_results -o fastqc_results

    echo "MultiQC report generated in fastqc_results/"
fi
