# Sequence_Alignment_Methods

# Computational Practice: Sequence Alignment on a High Performance Computer (HPC) cluster.

# Objective
In this lab, you will access a high-performance computing environment, install bioinformatics software, and perform sequence alignments. These are common computational biology tasks.

# Learning Objectives
* Access a high-performance computing environment
* Install and use different sequence alignment tools
* Compare global and local alignment algorithms
* Analyze alignment results for human genomic data

# Instructions
## Step 1: Access the Palmetto2 Cluster
Open your web browser and navigate to https://ondemand.rcd.clemson.edu/
Log in using your Clemson credentials
Complete any multi-factor authentication if required

## Step 2: Request Interactive Computing Resources
From the dashboard, click on "Interactive Apps" in the top navigation menu
Select "Jupyter Notebook" from the dropdown menu
Configure your session with the following specifications:
```
Number of cores: 4
Memory: 32GB
Wall time (job duration): 12:00:00 (12 hours)
Partition: Choose an appropriate partition based on availability
QOS: Select regular (or as recommended by your instructor)
```

Click "Launch" to submit your resource request
Wait for your job to start (may take a few minutes)

## Step 3: Access Terminal in Jupyter
When your Jupyter session is ready, click on "Connect to Jupyter"
In the Jupyter interface, click "New" in the top right corner
Select "Terminal" from the dropdown menu

## Step 4: Create Working Directory
```
# Navigate to your scratch space
cd /scratch/$USER

# Create a new directory for this project
mkdir -p sequence_alignment_project

# Navigate into the directory
cd sequence_alignment_project

# Verify your location
pwd
```

## Step 5: Install Required Software
```
# Load necessary anaconda modules (adjust based on Palmetto2's available modules). Anaconda environments let you install and isolate software in your home directory.  You can create multiple environments for different purposes.
module load anaconda3/2023.09-0 

# Create and activate  personal conda environment for sequence alignment tools
conda create -n alignment_env
source activate alignment_env

# Install BLAST (blastn,blastp,blastx,tblastn), Smith-Waterman (water) and Needleman-Wunsch (needle) sequence alignment tools
# EMBOSS package contains both water and needle software
conda install -c bioconda emboss

# Install BLAST+ suite
conda install -c bioconda blast
```
# Verify installations
```
# Check EMBOSS tools
water -help
needle -help
embossversion

# Check BLAST tools
blastn -help
blastp -help
blastx -help
tblastn -help
```

## Step 6: Download Human Genome and cDNA Datasets in your working directory on scratch
```
# Create directories for data in your working directory on scracth (not your home directory, please)
mkdir -p data/genome data/cdna

# Download human genome (GRCh38) chromosome 1 sequence 
cd /scratch/$USER/sequence_alignment_project/data/genome #Change $USER to your username
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#Note: The full hg38 genome is available from ENSEMBL at this URL: https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download human cDNA
cd /scratch/$USER/sequence_alignment_project/data/cdna
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Return to main project directory
cd /scratch/$USER/sequence_alignment_project/
```

## Step 7: Uncompress the Files
```
# Uncompress genome file
gunzip data/genome/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

# Uncompress cDNA file
gunzip data/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Verify uncompressed files
ls -lh data/genome/
ls -lh data/cdna/

# Look at the beginning of each file to confirm content
head -n 20 data/genome/Homo_sapiens.GRCh38.dna.chromosome.1.fa
head -n 20 data/cdna/Homo_sapiens.GRCh38.cdna.all.fa
```

## Step 8: Prepare for Alignment
```
# Extract a subset of the data for testing (full genome alignment would take too long)
# Extract first 100 cDNA sequences
head -n 1000 data/cdna/Homo_sapiens.GRCh38.cdna.chromosome.1.fa > data/cdna/test_cdna.fa

# Create BLAST database for the genome
cd data/genome/
makeblastdb -in Homo_sapiens.GRCh38.dna.chromosome.1.fa -dbtype nucl -out human_chr1
cd ../../
```

## Step 9: Perform Sequence Alignments
```
# Create directories for alignment results
mkdir -p results/blast results/smith_waterman

# Run BLAST alignment
blastn -query data/cdna/test_cdna.fa \
       -db data/genome/human_chr1 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -out results/blast/blast_alignment.tsv

# Run Smith-Waterman alignment using EMBOSS water
# Note: This will be much slower than BLAST, so we'll just do a few sequences
head -n 200 data/cdna/test_cdna.fa > data/cdna/mini_test.fa

# Run water (Smith-Waterman) on a single cDNA sequence against a small portion of the genome
head -n 2 data/cdna/mini_test.fa > data/cdna/single_cdna.fa
head -n 10000 data/genome/Homo_sapiens.GRCh38.dna.chromosome.1.fa > data/genome/mini_genome.fa

water -asequence data/cdna/single_cdna.fa \
      -bsequence data/genome/mini_genome.fa \
      -gapopen 10 -gapextend 0.5 \
      -outfile results/smith_waterman/water_alignment.txt
```
## Step 10: Compare Alignment Results
```
# Look at BLAST results
head results/blast/blast_alignment.tsv

# Count the number of alignments found by BLAST
wc -l results/blast/blast_alignment.tsv

# Check Smith-Waterman alignment
cat results/smith_waterman/water_alignment.txt

# Create a Python script to analyze the results
cat << 'EOF' > compare_alignments.py
#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import os

# Read BLAST results
blast_results = pd.read_csv("results/blast/blast_alignment.tsv", sep='\t', 
                           names=["qseqid", "sseqid", "pident", "length", "mismatch", 
                                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

# Basic statistics on BLAST results
print("BLAST Alignment Statistics:")
print(f"Total alignments: {len(blast_results)}")
print(f"Average percent identity: {blast_results['pident'].mean():.2f}%")
print(f"Average alignment length: {blast_results['length'].mean():.2f}")
print(f"Average E-value: {blast_results['evalue'].mean():.2e}")

# Create results directory for plots
os.makedirs("results/plots", exist_ok=True)

# Create histogram of percent identities
plt.figure(figsize=(10, 6))
plt.hist(blast_results['pident'], bins=20, alpha=0.7)
plt.title("Distribution of Alignment Percent Identities")
plt.xlabel("Percent Identity")
plt.ylabel("Count")
plt.savefig("results/plots/percent_identity_dist.png")

# Create scatter plot of alignment length vs. percent identity
plt.figure(figsize=(10, 6))
plt.scatter(blast_results['length'], blast_results['pident'], alpha=0.5)
plt.title("Alignment Length vs. Percent Identity")
plt.xlabel("Alignment Length")
plt.ylabel("Percent Identity")
plt.savefig("results/plots/length_vs_identity.png")

# Analysis of Smith-Waterman results would typically compare specific alignments
# but since we only ran it on a small subset, we'll just note how it would be done
print("\nSmith-Waterman vs BLAST comparison:")
print("For a complete comparison, you would need to:")
print("1. Run both algorithms on the same sequences")
print("2. Compare alignment positions, gaps, and scores")
print("3. Evaluate sensitivity and specificity")
print("4. Compare computational performance (time and memory usage)")

EOF

# Make the script executable
chmod +x compare_alignments.py

# Run the analysis
python compare_alignments.py
```

# Questions for Reflection
* What are the key differences you observed between BLAST and Smith-Waterman alignments?
* Why is BLAST significantly faster than Smith-Waterman for genomic sequence alignments?
* In what scenarios would you prefer to use Smith-Waterman over BLAST?
* How does the E-value in BLAST relate to the statistical significance of an alignment?
* What challenges did you encounter aligning cDNA sequences to genomic DNA, and why?
* How would the alignment results differ if you used Needleman-Wunsch instead of Smith-Waterman?
* Based on your results, what percentage of your cDNA sequences aligned well to the genome?

# Generative AI Prompts
## Prompt 1: Understanding Pairwise vs. Global Sequence Alignment
```
I'm working with sequence alignment algorithms and need to understand the conceptual and practical differences between pairwise and global alignment approaches. Could you explain:

1. The fundamental differences between local, global, and semi-global sequence alignment
2. When each approach is most appropriate in bioinformatics workflows
3. How scoring matrices and gap penalties affect each alignment type differently
4. The mathematical foundations of these alignment strategies (dynamic programming approaches)
5. How these approaches handle evolutionary events like insertions, deletions, and substitutions
6. Real-world examples where choosing the wrong alignment strategy would lead to incorrect biological conclusions
7. How modern alignment tools implement these different strategies
8. Computational complexity comparisons between different alignment types
9. Visualization techniques for comparing alignment results from different methods
10. How these fundamental alignment concepts extend to multiple sequence alignment
```

## Prompt 2: Comparing BLAST and Smith-Waterman Algorithms
```
I'm comparing sequence alignment results between BLAST and Smith-Waterman for a genomics project. Could you provide a detailed explanation of:

1. The core algorithmic differences between BLAST and Smith-Waterman
2. Why BLAST is faster but Smith-Waterman is more sensitive (with specific examples)
3. The heuristics BLAST uses and what trade-offs they introduce
4. How BLAST's seed-and-extend approach differs from Smith-Waterman's dynamic programming matrix
5. The statistical significance measures used in each algorithm (E-values vs. raw scores)
6. Scenarios where BLAST might miss alignments that Smith-Waterman would detect
7. How parameter choices (match/mismatch scores, gap penalties) affect each algorithm differently
8. Modern optimizations of Smith-Waterman (GPU implementations, SIMD instructions) and how they affect performance
9. The impact of sequence length and composition on the performance of each algorithm
10. Best practices for comparing and integrating results from both algorithms for maximum biological insight
```

# Additional Resources
* BLAST Documentation (https://www.ncbi.nlm.nih.gov/books/NBK279690/)
* EMBOSS Water Documentation (https://www.bioinformatics.nl/cgi-bin/emboss/help/water)
* EMBOSS Needle Documentation (https://www.bioinformatics.nl/cgi-bin/emboss/help/needle)
* Palmetto2 User Guide (https://docs.rcd.clemson.edu/palmetto/index.html)
