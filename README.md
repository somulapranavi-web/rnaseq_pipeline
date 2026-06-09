# Simulated RNA-Seq Analysis Pipeline

An end-to-end RNA-Seq bioinformatics analysis pipeline. It simulates synthetic paired-end sequencing reads across a 3 vs. 3 (Control vs. Treated) experimental design, performs quality control, quantifies transcript abundances, and conducts differential gene expression analysis.

## Pipeline Architecture

```text
               +-----------------------------+
               |  generate_data.py (Python)  |
               +--------------+--------------+
                              | (Simulated FASTQs)
                              v
                 +------------+------------+
                 |   FastQC & MultiQC (QC) |
                 +------------+------------+
                              | (Quality verification)
                              v
                  +-----------+-----------+
                  |   Salmon Quant (Align) |
                  +-----------+-----------+
                              | (Abundance values)
                              v
                +-------------+-------------+
                |     DESeq2 Analysis (R)   |
                +-------------+-------------+
                              |
                     +--------+--------+
                     |                 |
                     v                 v
            diff_expr_results.csv   volcano_plot.png
```

---

## File Structure

- `data/`: Location for raw input FASTQ read files (ignored by git; generated via script).
- `reference/`: Contains `transcripts.fasta` (reference transcriptome).
- `scripts/`:
  - `generate_data.py`: Generates the mock reference and mock read datasets.
  - `setup_env.sh`: Automatically configures the WSL environment and Conda packages.
  - `run_pipeline.sh`: Runs QC (FastQC & MultiQC) and quantification (Salmon).
  - `diff_expr.R`: Conducts statistical analysis (DESeq2) and plots the results.
- `diff_expr_results.csv`: Table of gene-wise log2 fold changes and p-values.
- `volcano_plot.png`: Graphical output highlighting statistically significant genes.

---

## Getting Started

The pipeline is configured to run inside a **WSL (Windows Subsystem for Linux)** Ubuntu environment.

### 1. Generate Synthetic Data
First, run the generator script using Python 3 to simulate 6 paired-end samples (Control1-3, Treated1-3, 25k read pairs each):
```bash
python scripts/generate_data.py
```

### 2. Setup the WSL Environment
Run the setup script inside WSL to install Conda, configure channels, and build the `rnaseq_env` environment:
```bash
bash scripts/setup_env.sh
```

### 3. Run the Pipeline
Execute the main pipeline script to run QC and Salmon quantification:
```bash
bash scripts/run_pipeline.sh
```

### 4. Differential Gene Expression Analysis
Execute the R script to run DESeq2, write out the results CSV, and plot the Volcano Plot:
```bash
# Make sure the rnaseq_env is active
source ~/miniconda3/etc/profile.d/conda.sh && conda activate rnaseq_env
Rscript scripts/diff_expr.R
```

---

## Simulated Biology Details

During simulation (`generate_data.py`), we target specific genes for differential expression:
- **Upregulated Targets (5x change):** `YAL005W`, `YAL015W`, `YAL025W`
- **Downregulated Targets (0.1x change):** `YAL010W`, `YAL020W`, `YAL030W`

The DESeq2 results accurately identify these targets with statistically significant adjusted p-values ($P_{\text{adj}} \approx 0$).
