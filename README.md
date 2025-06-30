# DNA_seq_worm_microbiome_taxonomy_classification

# Necator americanus Microbiome Analysis Pipeline

**This is your README.md file - copy this entire content into a file called README.md in your repository**

A comprehensive bioinformatics pipeline for analyzing bacterial communities associated with *Necator americanus* (hookworm) samples using unmapped sequencing reads.

## Overview

This pipeline extracts unmapped reads from aligned BAM files, performs taxonomic classification using Kraken2/Bracken, and generates publication-ready visualizations of bacterial community composition. The workflow is optimized for SLURM cluster environments.

## Pipeline Workflow

```
BAM files → Unmapped reads → QC → Trimming → Kraken2 → Bracken → Filtering → Visualization
```

### Steps:
1. **Extract unmapped reads** from aligned BAM files
2. **Quality control** with FastQC/MultiQC
3. **Read trimming** with Trimmomatic
4. **Taxonomic classification** with Kraken2
5. **Abundance estimation** with Bracken
6. **Bacterial filtering** and normalization
7. **Statistical analysis** and visualization

## Requirements

### Software Dependencies
- **SLURM** workload manager
- **samtools** (≥1.10)
- **FastQC** (≥0.11.9)
- **MultiQC** (≥1.9)
- **Trimmomatic** (≥0.39)
- **Kraken2** (v2.0.8-beta, for Bracken compatibility)
- **Bracken** (≥2.6.0)
- **Python** (≥3.7) with packages:
  - pandas
  - numpy
  - matplotlib
  - collections
- **R** (≥4.0) with packages:
  - vegan
  - ecodist

### Database Requirements
- **Kraken2 database** (set `$KRAKEN2_DB` environment variable)
  - Recommended: Standard database or custom database with bacterial genomes

## Installation

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/necator-microbiome-pipeline.git
cd necator-microbiome-pipeline
```

2. **Set up environment:**
```bash
# Load required modules (adjust for your HPC system)
module load samtools fastqc multiqc trimmomatic kraken2/2.0.8-beta bracken

# Set Kraken2 database path
export KRAKEN2_DB=/path/to/your/kraken2/database
```

3. **Install R packages:**
```bash
Rscript install_packages.R
```

## Input Data

- **BAM files**: Aligned sequencing reads (`.sorted.refrename.bam`)
- Expected naming convention: `sample_name.sorted.refrename.bam`
- Files should be placed in the working directory

## Usage

### Quick Start
```bash
# 1. Extract unmapped reads (adjust array size for your sample count)
sbatch --array=1-N 1_run_get_unmapped.sh

# 2. Quality control
sbatch --array=1-18 2x_run_fastqc.sh  # Adjust array size for fastq count

# 3. Trim reads
sbatch --array=1-N 3x_run_trimming.sh

# 4. Taxonomic classification
sbatch --array=1-N 4x_run_kraken2.sh

# 5. Combine and filter bacterial data
sbatch 5x_combine_filter_bacteria.sh

# 6. Normalize and create pie chart data
sbatch 6x_normalize_combine_bacteria.sh

# 7. Create visualizations
python3 make_pie_chart.py

# 8. Statistical analysis (optional)
Rscript microbiome_analysis.R
```

### Detailed Steps

#### Step 1: Extract Unmapped Reads
```bash
sbatch --array=1-9 1_run_get_unmapped.sh
```
- Extracts paired unmapped reads from BAM files
- Outputs: `fastq_for_assembly/sample.unmapped.R1.fastq` and `R2.fastq`

#### Step 2: Quality Control
```bash
sbatch --array=1-18 2x_run_fastqc.sh
```
- Runs FastQC on all FASTQ files
- Generates MultiQC summary report
- Output: `fastqc_results/` directory

#### Step 3: Read Trimming
```bash
sbatch --array=1-9 3x_run_trimming.sh
```
- Trims low-quality bases and short reads
- Output: `trimmed_fastq/` directory

#### Step 4: Taxonomic Classification
```bash
sbatch --array=1-9 4x_run_kraken2.sh
```
- Classifies reads with Kraken2
- Estimates abundances with Bracken at multiple taxonomic levels
- Outputs: `kraken2_output/` and `bracken_output/` directories

#### Step 5: Bacterial Data Processing
```bash
sbatch 5x_combine_filter_bacteria.sh
```
- Extracts bacterial genera from all samples
- Combines into count matrices
- Outputs: 
  - `bacterial_genera_raw_counts.txt`
  - `bacterial_genera_relative_abundance.txt`
  - `bacterial_matrix_summary.txt`

#### Step 6: Normalization and Filtering
```bash
sbatch 6x_normalize_combine_bacteria.sh
```
- Performs rarefaction normalization
- Applies statistical filtering
- Creates pie chart-ready data
- Outputs:
  - `bacteria_text_files/bacterial_genera_rarefied_counts_filtered.txt`
  - `bacteria_text_files/bacterial_overall_composition_for_piechart_new.txt`

#### Step 7: Visualization
```bash
python3 make_pie_chart.py
```
- Creates publication-ready pie charts and bar charts
- Outputs: PDF, PNG, and SVG formats

## Output Files

### Primary Results
- **`bacterial_community_pie_chart.pdf`** - Publication-ready pie chart
- **`bacterial_genera_rarefied_counts_filtered.txt`** - Normalized count matrix
- **`bacterial_normalization_summary_new.txt`** - Analysis summary

### Intermediate Files
- **`fastqc_results/`** - Quality control reports
- **`trimmed_fastq/`** - Quality-trimmed reads
- **`kraken2_output/`** - Taxonomic classification results
- **`bracken_output/`** - Abundance estimation results
- **`bacteria_text_files/`** - Processed bacterial data

### Statistical Analysis (Optional)
- **`NMDS_plot.pdf`** - NMDS ordination plot
- **`diversity_indices.csv`** - Alpha diversity metrics
- **`community_matrix.csv`** - Species-level community matrix

## Configuration

### Adjusting for Your Data

1. **Sample count**: Modify array parameters in SLURM scripts
   ```bash
   #SBATCH --array=1-N  # Change N to your sample count
   ```

2. **Resource allocation**: Adjust memory and time limits based on data size
   ```bash
   #SBATCH --mem=50G     # Increase for larger datasets
   #SBATCH --time=4:00:00
   ```

3. **Filtering parameters**: Edit filtering thresholds in `normalize_bacterial_abundance.py`
   ```python
   min_reads = 2000         # Minimum read count threshold
   min_abundance_pct = 1.0  # Minimum abundance percentage
   ```

## Quality Control and Validation

### Data Quality Checks
- **FastQC/MultiQC reports**: Check read quality metrics
- **Kraken2 classification rates**: Typically 10-30% for microbiome samples
- **Rarefaction depth**: Ensures fair comparison across samples
- **Filtering statistics**: Removes low-confidence assignments

### Expected Results
- **Classification rate**: 15-40% of unmapped reads classified
- **Bacterial diversity**: 20-100 genera depending on sample type
- **Dominant genera**: Should reflect known microbiome patterns

## Troubleshooting

### Common Issues

1. **Low classification rates**:
   - Check Kraken2 database completeness
   - Verify input data quality
   - Consider host contamination removal

2. **Memory errors**:
   - Increase memory allocation in SLURM scripts
   - Use smaller Kraken2 database if available

3. **No bacterial reads detected**:
   - Verify samples contain microbial DNA
   - Check host genome alignment completeness
   - Review Kraken2 database bacterial content

### File Paths
Update absolute paths in scripts to match your system:
```bash
cd /path/to/your/working/directory  # Update in all scripts
```

## Citation

If you use this pipeline, please cite:
- **Kraken2**: Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019).
- **Bracken**: Lu, J., Breitwieser, F.P., Thielen, P. & Salzberg, S.L. Bracken: estimating species abundance in metagenomics data. PeerJ Comput. Sci. 3:e104 (2017).

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Contact

- **Author**: Zach Pella
- **Institution**: Fauver Lab
- **Issues**: Please use GitHub issues for bug reports and feature requests

---

## File Structure
```
necator-microbiome-pipeline/
├── README.md
├── LICENSE
├── install_packages.R
├── scripts/
│   ├── 1_run_get_unmapped.sh
│   ├── 2x_run_fastqc.sh
│   ├── 3x_run_trimming.sh
│   ├── 4x_run_kraken2.sh
│   ├── 5x_combine_filter_bacteria.sh
│   ├── 6x_normalize_combine_bacteria.sh
│   ├── make_pie_chart.py
│   ├── normalize_bacterial_abundance.py
│   └── microbiome_analysis.R
├── config/
│   └── pipeline_config.yaml
├── docs/
│   ├── installation.md
│   ├── tutorial.md
│   └── troubleshooting.md
└── example_data/
    └── README.md
```
