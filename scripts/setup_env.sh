#!/bin/bash
# setup_env.sh
# Automates the installation of Miniconda and NGS tools inside WSL Ubuntu.

set -e # Exit immediately on erro

echo "=== Starting RNA-Seq Environment Setup ==="

# Define Miniconda path
MINICONDA_DIR="$HOME/miniconda3"
CONDA_SH="$MINICONDA_DIR/etc/profile.d/conda.sh"

# 1. Install Miniconda if not present
if [ ! -f "$CONDA_SH" ]; then
    echo "Miniconda not found. Downloading and installing..."
    curl -L -o miniconda_installer.sh "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    bash miniconda_installer.sh -b -p "$MINICONDA_DIR"
    rm miniconda_installer.sh
    echo "Miniconda installed successfully."
else
    echo "Miniconda is already installed at $MINICONDA_DIR."
fi

# 2. Source conda configuration
source "$CONDA_SH"
conda init bash --user || true

# 3. Configure channels for Bioinformatics (Conda Forge and Bioconda)
echo "Configuring Conda channels..."
conda config --remove channels defaults || true
conda config --add channels conda-forge || true
conda config --add channels bioconda || true
conda config --set channel_priority strict || true
conda config --set solver libmamba || true

# 4. Create conda environment if it doesn't exist
ENV_NAME="rnaseq_env"
ENV_PATH="$MINICONDA_DIR/envs/$ENV_NAME"

if [ -d "$ENV_PATH" ]; then
    echo "Conda environment '$ENV_NAME' already exists at $ENV_PATH."
else
    echo "Creating environment '$ENV_NAME' using micromamba solver..."
    
    # Download micromamba binary using python standard library to avoid bzip2 dependencies
    mkdir -p "$HOME/bin"
    python3 -c "import urllib.request, tarfile, io, os; url = 'https://micro.mamba.pm/api/micromamba/linux-64/latest'; req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla'}); response = urllib.request.urlopen(req); tf = tarfile.open(fileobj=io.BytesIO(response.read()), mode='r:bz2'); tf.extract('bin/micromamba', path=os.path.expanduser('~'))"
    
    # Run micromamba to create the environment directly in miniconda envs directory
    "$HOME/bin/micromamba" create -p "$ENV_PATH" -c conda-forge -c bioconda --channel-priority flexible -y \
        salmon \
        fastqc \
        multiqc \
        r-base \
        r-ggplot2 \
        r-matrix \
        bioconductor-deseq2
        
    echo "Conda environment '$ENV_NAME' created successfully."
fi

echo "=== Setup Completed Successfully! ==="
echo "You can activate this environment inside WSL using:"
echo "  source $CONDA_SH && conda activate $ENV_NAME"
