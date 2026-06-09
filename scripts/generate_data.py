#!/usr/bin/env python3
"""
Simulated RNA-Seq Data Generator
Generates mock reference transcripts and paired-end FastQ files for a 3 vs 3 design.
"""

import os
import random
import sys

# Set random seed for reproducibility
random.seed(42)

# Vocabulary of nucleotides
NUCLEOTIDES = ['A', 'C', 'G', 'T']

def generate_random_sequence(length: int) -> str:
    """Generates a random DNA sequence of specified length."""
    return "".join(random.choices(NUCLEOTIDES, k=length))

def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(comp.get(base, 'N') for base in reversed(seq))

def main():
    print("=== Generating Simulated RNA-Seq Dataset ===")
    
    # Establish directories
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    ref_dir = os.path.join(base_dir, "reference")
    
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(ref_dir, exist_ok=True)
    
    # 1. Define transcripts (50 genes of lengths between 500 and 1500 bp)
    num_genes = 50
    genes = {}
    gene_lengths = {}
    for i in range(1, num_genes + 1):
        gene_id = f"YAL{i:03d}W"  # Yeast systematic naming format
        length = random.randint(500, 1500)
        genes[gene_id] = generate_random_sequence(length)
        gene_lengths[gene_id] = length
        
    # Write reference transcripts to FASTA
    ref_fasta_path = os.path.join(ref_dir, "transcripts.fasta")
    with open(ref_fasta_path, 'w') as f:
        for gene_id, seq in genes.items():
            f.write(f">{gene_id}\n{seq}\n")
    print(f"Created reference transcripts: {ref_fasta_path} ({num_genes} genes)")
    
    # 2. Define expression profiles
    # Control abundance baseline (exponential decay to simulate biological variation)
    base_abundances = {gene_id: (0.9 ** idx) for idx, gene_id in enumerate(genes.keys())}
    # Normalize abundances to sum to 1.0
    total_base = sum(base_abundances.values())
    base_abundances = {k: v / total_base for k, v in base_abundances.items()}
    
    # Define differential expression target genes
    # 3 upregulated genes, 3 downregulated genes
    upregulated_targets = ["YAL005W", "YAL015W", "YAL025W"]
    downregulated_targets = ["YAL010W", "YAL020W", "YAL030W"]
    
    # Generate 6 samples (3 Control vs 3 Treated)
    samples = {
        "Control1": "control",
        "Control2": "control",
        "Control3": "control",
        "Treated1": "treated",
        "Treated2": "treated",
        "Treated3": "treated"
    }
    
    reads_per_sample = 25000  # Number of read-pairs to generate (50k total reads per sample)
    read_length = 75
    fragment_mean = 250
    fragment_std = 30
    
    for sample_name, condition in samples.items():
        print(f"Simulating reads for {sample_name} ({condition})...")
        
        # Adjust abundances based on condition
        sample_abundances = base_abundances.copy()
        
        if condition == "treated":
            # Upregulate targets by 5x
            for target in upregulated_targets:
                sample_abundances[target] *= 5.0
            # Downregulate targets by 0.1x (90% reduction)
            for target in downregulated_targets:
                sample_abundances[target] *= 0.1
                
            # Renormalize
            tot = sum(sample_abundances.values())
            sample_abundances = {k: v / tot for k, v in sample_abundances.items()}
            
        # Cumulative probabilities for sampling
        gene_list = list(sample_abundances.keys())
        weights = list(sample_abundances.values())
        
        # Files to write
        r1_path = os.path.join(data_dir, f"{sample_name}_1.fastq")
        r2_path = os.path.join(data_dir, f"{sample_name}_2.fastq")
        
        # Phred score (F = Phred 37, standard high quality)
        qual_str = "F" * read_length
        
        with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
            # We sample reads with replacement
            sampled_genes = random.choices(gene_list, weights=weights, k=reads_per_sample)
            
            for idx, gene_id in enumerate(sampled_genes):
                gene_seq = genes[gene_id]
                gene_len = gene_lengths[gene_id]
                
                # Determine fragment size
                frag_len = int(random.normalvariate(fragment_mean, fragment_std))
                frag_len = max(read_length + 20, min(frag_len, gene_len))
                
                # Determine fragment start position
                max_start = gene_len - frag_len
                start = random.randint(0, max_start)
                
                # Extract fragment
                fragment = gene_seq[start:start + frag_len]
                
                # Generate Read 1 (5' end of fragment)
                read1 = fragment[:read_length]
                # Introduce occasional sequencing errors (0.1% rate)
                read1_list = list(read1)
                for pos in range(read_length):
                    if random.random() < 0.001:
                        read1_list[pos] = random.choice([n for n in NUCLEOTIDES if n != read1_list[pos]])
                read1 = "".join(read1_list)
                
                # Generate Read 2 (3' end of fragment, reverse-complemented)
                read2_frag = fragment[-read_length:]
                read2 = reverse_complement(read2_frag)
                read2_list = list(read2)
                for pos in range(read_length):
                    if random.random() < 0.001:
                        read2_list[pos] = random.choice([n for n in NUCLEOTIDES if n != read2_list[pos]])
                read2 = "".join(read2_list)
                
                # Write Read 1 in FASTQ format
                r1.write(f"@read_{idx}_{gene_id}\n")
                r1.write(f"{read1}\n")
                r1.write("+\n")
                r1.write(f"{qual_str}\n")
                
                # Write Read 2 in FASTQ format
                r2.write(f"@read_{idx}_{gene_id}\n")
                r2.write(f"{read2}\n")
                r2.write("+\n")
                r2.write(f"{qual_str}\n")
                
        print(f"  Saved FASTQ files: {sample_name}_1.fastq / {sample_name}_2.fastq")
        
    print("=== Data Generation Completed Successfully ===")

if __name__ == "__main__":
    main()
