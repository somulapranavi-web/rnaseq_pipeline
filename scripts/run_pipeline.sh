#!/bin/bash
# run_pipeline.sh
# Runs the RNA-Seq pipeline: FastQC -> MultiQC -> Salmon Index -> Salmon Quant.

set -e # Exit immediately on error

echo "=== Starting RNA-Seq Analysis Pipeline ==="

# Define Miniconda paths
MINICONDA_DIR="$HOME/miniconda3"
CONDA_SH="$MINICONDA_DIR/etc/profile.d/conda.sh"

if [ -f "$CONDA_SH" ]; then
    source "$CONDA_SH"
else
    echo "Error: Miniconda not found at $MINICONDA_DIR. Please run setup_env.sh first."
    exit 1
fi

# Activate the conda environment
ENV_NAME="rnaseq_env"
echo "Activating Conda environment '$ENV_NAME'..."
conda activate "$ENV_NAME"

# Establish directories
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPTS_DIR/.." && pwd)"
DATA_DIR="$BASE_DIR/data"
REF_DIR="$BASE_DIR/reference"
QC_DIR="$BASE_DIR/qc"
QUANTS_DIR="$BASE_DIR/quants"

mkdir -p "$QC_DIR"
mkdir -p "$QUANTS_DIR"

# 1. Run Quality Control (FastQC)
echo "=== Step 1: Running FastQC Quality Control ==="
fastqc -o "$QC_DIR" "$DATA_DIR"/*.fastq

# 2. Run MultiQC to aggregate report
echo "=== Step 2: Running MultiQC compilation ==="
multiqc -o "$QC_DIR" --force "$QC_DIR"

# 3. Build Salmon index
echo "=== Step 3: Building Salmon index ==="
SALMON_INDEX="$REF_DIR/salmon_index"
if [ -d "$SALMON_INDEX" ]; then
    echo "Salmon index already exists. Skipping indexing."
else
    salmon index -t "$REF_DIR/transcripts.fasta" -i "$SALMON_INDEX" -k 31
    echo "Salmon indexing complete."
fi

# 4. Run Salmon quantification
echo "=== Step 4: Quantifying transcript abundances ==="
for sample in Control1 Control2 Control3 Treated1 Treated2 Treated3; do
    if [ -f "$DATA_DIR/${sample}_1.fastq" ] && [ -f "$DATA_DIR/${sample}_2.fastq" ]; then
        echo "Quantifying sample: $sample..."
        salmon quant -i "$SALMON_INDEX" -l A \
            -1 "$DATA_DIR/${sample}_1.fastq" \
            -2 "$DATA_DIR/${sample}_2.fastq" \
            -p 4 \
            --validateMappings \
            -o "$QUANTS_DIR/${sample}_quant"
    else
        echo "Warning: FastQ files for $sample not found. Skipping."
    fi
done

echo "=== Pipeline Completed Successfully! ==="
echo "Quality control results: $QC_DIR/multiqc_report.html"
echo "Quantification results: $QUANTS_DIR/"
