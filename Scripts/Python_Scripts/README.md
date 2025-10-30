# MK Table Builder 

## Purpose
Prepares data for McDonald-Kreitman (MK) test or Direction of Selection (DoS) calculation by classifying SNPs in coding sequences and building contingency tables for detecting selection.

## Installation

### Create conda environment from YAML file:
```bash
# Create the environment
conda env create -f mk_table_builder_env.yml

# Activate the environment
conda activate mk_table_builder
```

## What it does
1. **Extracts coding sequences (CDS)** from a genome using GFF annotations
2. **Identifies SNPs** in coding regions from a VCF file
3. **Classifies each SNP** as:
   - **Synonymous (S)** or **Nonsynonymous/Replacement (R)** - does it change the amino acid?
   - **Polymorphic (P)** or **Fixed (F)** - is it variable within the ingroup or fixed between ingroup/outgroup?
4. **Counts SNPs per gene** to create a 2×2 contingency table:
   ```
                Polymorphic | Fixed
   Synonymous      PS      |  FS
   Nonsynonymous   PR      |  FR
   ```
5. **Calculates Tsil and Trepl** - expected synonymous and nonsynonymous sites based on codon usage

## Output
CSV file with columns: `geneID, PR, FR, PS, FS, Tsil, Trepl, nout, npop`
- **PR**: Polymorphic Nonsynonymous (Replacement)
- **FR**: Fixed Nonsynonymous (Replacement)  
- **PS**: Polymorphic Synonymous
- **FS**: Fixed Synonymous
- **Tsil**: Expected synonymous sites
- **Trepl**: Expected replacement sites
- **nout**: Number of outgroup chromosomes
- **npop**: Number of ingroup chromosomes


## Input requirements
- **Genome FASTA**: Reference genome sequences
- **GFF3**: Gene annotations with CDS features
- **VCF**: Variant calls with genotypes for all samples
- **Ingroup samples**: Text file listing ingroup sample IDs (one per line)
- **Outgroup samples**: Text file listing outgroup sample IDs (one per line)

## Usage
```bash
# Activate the conda environment
conda activate mk_table_builder

# Run the script
python Build-mkTest-Table.py \
  -s genome.fasta \
  -g annotations.gff3 \
  -v variants.vcf \
  -ig ingroup_samples.txt \
  -og outgroup_samples.txt \
  -o mk_results.csv
```
