def ensure_vcf_index(vcf_file):
    """
    Check if VCF index file exists, and create it if it doesn't.
    Returns True if indexing was successful or index already exists.
    """
    # Check for .tbi (tabix) index for gzipped VCF
    index_file = vcf_file + ".tbi"
    if os.path.exists(index_file):
        print(f"VCF index file {index_file} already exists.")
        return True
    
    # Check for .csi index (default for bcftools)
    csi_index = vcf_file + ".csi"
    if os.path.exists(csi_index):
        print(f"VCF index file {csi_index} already exists.")
        return True
    
    # Create index
    print(f"Creating index for {vcf_file}...")
    
    # Determine if file is bgzipped by checking file signature
    is_bgzipped = False
    try:
        with open(vcf_file, 'rb') as f:
            # Check for gzip magic number (first two bytes)
            if f.read(2) == b'\x1f\x8b':
                is_bgzipped = True
    except Exception as e:
        print(f"Error checking if file is gzipped: {e}")
        return False
    
    # Index command depends on whether file is compressed
    if is_bgzipped:
        cmd = f"bcftools index {vcf_file}"
    else:
        # For uncompressed VCF, we need to compress it first
        print("VCF file is not compressed. Compressing first...")
        compressed_vcf = vcf_file + ".gz"
        compress_cmd = f"bgzip -c {vcf_file} > {compressed_vcf}"
        result = run_command(compress_cmd)
        if not result:
            print("Compression failed.")
            return False
        
        # Update vcf_file to point to the compressed file
        vcf_file = compressed_vcf
        cmd = f"bcftools index {vcf_file}"
    
    # Run the indexing command
    result = run_command(cmd)
    if result is None:
        print("Indexing failed.")
        return False
    
    print("Indexing completed successfully.")
    return True#!/usr/bin/env python3
"""
Extract SNPs from specific genes in a VCF file and save to text files.
This script extracts variants from four genes and saves only the SNPs to separate text files.
It uses gene FASTA files to accurately predict amino acid changes in the format XposY.
"""

import subprocess
import os
import sys
import re
from collections import defaultdict

# Genetic code - DNA codon to amino acid mapping
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def run_command(command):
    """Run a shell command and return the output."""
    print(f"Running: {command}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command: {result.stderr}")
        return None
    return result.stdout

def read_fasta(fasta_file):
    """Read a FASTA file and return the sequence."""
    if not os.path.exists(fasta_file):
        print(f"Warning: FASTA file {fasta_file} not found.")
        return None
        
    sequence = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence += line.strip()
    return sequence.upper()

def translate_dna(dna_sequence):
    """Translate a DNA sequence to protein."""
    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:  # Make sure we have a complete codon
            amino_acid = GENETIC_CODE.get(codon, 'X')  # 'X' for unknown codons
            protein += amino_acid
    return protein

def get_amino_acid_change(gene_sequence, position, ref, alt, gene_start):
    """
    Calculate the amino acid change based on the genomic position and the gene sequence.
    
    Args:
        gene_sequence: The coding sequence of the gene
        position: Genomic position of the variant
        ref: Reference allele
        alt: Alternate allele
        gene_start: Start position of the gene in the genome
        
    Returns:
        String representing amino acid change in XposY format, or None if not a coding change
    """
    # Calculate position within the gene sequence
    relative_position = position - gene_start
    
    # Check if position is within the gene sequence
    if relative_position < 0 or relative_position >= len(gene_sequence):
        return "non-coding"
    
    # Find the codon position and offset
    codon_position = relative_position // 3
    codon_offset = relative_position % 3
    
    # Extract the affected codon
    start_idx = codon_position * 3
    if start_idx + 2 >= len(gene_sequence):
        return "incomplete-codon"
        
    ref_codon = gene_sequence[start_idx:start_idx+3]
    
    # Create the variant codon
    alt_codon = list(ref_codon)
    alt_codon[codon_offset] = alt
    alt_codon = ''.join(alt_codon)
    
    # Translate codons to amino acids
    ref_aa = GENETIC_CODE.get(ref_codon, 'X')
    alt_aa = GENETIC_CODE.get(alt_codon, 'X')
    
    # Return the change (even for synonymous mutations)
    if ref_aa != alt_aa:
        return f"{ref_aa}{codon_position+1}{alt_aa}"
    else:
        return f"{ref_aa}{codon_position+1}{alt_aa} (synonymous)"

def load_gene_info():
    """
    Load gene information from FASTA files.
    Returns a dictionary with gene details.
    """
    genes = {
        "NPSR1": {
            "fasta": "assets/npsr1.fasta",
            "hg19": {"chr": "7", "start": 34694197, "end": 34885439},
            "hg38": {"chr": "7", "start": 34658413, "end": 34849655},
            "strand": "+"
        },
        "ADRB1": {
            "fasta": "assets/adrb1.fasta",
            "hg19": {"chr": "10", "start": 115804309, "end": 115805742},
            "hg38": {"chr": "10", "start": 114044133, "end": 114045566},
            "strand": "+"
        },
        "DEC2": {
            "fasta": "assets/dec2.fasta",
            "hg19": {"chr": "12", "start": 26125956, "end": 26128669},
            "hg38": {"chr": "12", "start": 26122066, "end": 26124779},
            "strand": "-"
        },
        "GRM1": {
            "fasta": "assets/grm1.fasta",
            "hg19": {"chr": "6", "start": 146755366, "end": 147160172},
            "hg38": {"chr": "6", "start": 146029518, "end": 146434796},
            "strand": "+"
        }
    }
    
    # Load sequences from FASTA files
    for gene_name, gene_data in genes.items():
        fasta_file = gene_data["fasta"]
        sequence = read_fasta(fasta_file)
        if sequence:
            gene_data["sequence"] = sequence
        else:
            print(f"Warning: Could not load sequence for {gene_name}")
            gene_data["sequence"] = None
    
    return genes

def extract_snps_with_aa_changes(vcf_file, use_hg38=True, chr_prefix=""):
    """
    Extract SNPs for each gene and annotate with amino acid changes.
    
    Args:
        vcf_file: Path to the VCF file
        use_hg38: Whether to use hg38 coordinates (True) or hg19 (False)
        chr_prefix: Chromosome prefix ("chr" or "")
    """
    # Create output directory
    output_dir = "gene_variants"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load gene information
    genes = load_gene_info()
    
    # Process each gene
    for gene_name, gene_data in genes.items():
        # Get genome build specific coordinates
        build = "hg38" if use_hg38 else "hg19"
        chrom = chr_prefix + gene_data[build]["chr"]
        start = gene_data[build]["start"]
        end = gene_data[build]["end"]
        sequence = gene_data.get("sequence")
        strand = gene_data.get("strand", "+")
        
        print(f"\nProcessing {gene_name} at {chrom}:{start}-{end}")
        
        if not sequence:
            print(f"Skipping {gene_name} - sequence not available")
            continue
            
        # Extract the gene's SNPs to a temporary file
        temp_vcf = os.path.join(output_dir, f"{gene_name}_temp.vcf")
        
        # Use indexed access if possible (faster)
        is_indexed = os.path.exists(vcf_file + ".tbi") or os.path.exists(vcf_file + ".csi")
        if is_indexed:
            extract_cmd = f"bcftools view -r {chrom}:{start}-{end} {vcf_file} > {temp_vcf}"
        else:
            # Fallback to slower grep-based extraction
            extract_cmd = f"bcftools view {vcf_file} | grep -P '^{chrom}\\t[0-9]+\\t' | awk '{{if ($2 >= {start} && $2 <= {end}) print}}' > {temp_vcf}"
        
        run_command(extract_cmd)
        
        if not os.path.exists(temp_vcf) or os.path.getsize(temp_vcf) == 0:
            print(f"No variants found for {gene_name} or extraction failed.")
            continue
            
        # For negative strand genes, we need to reverse complement the sequence
        working_sequence = sequence
        if strand == "-":
            working_sequence = reverse_complement(sequence)
        
        # Process SNPs and add amino acid changes
        output_file = os.path.join(output_dir, f"{gene_name}_snps.txt")
        snp_count = 0
        
        with open(temp_vcf, 'r') as vcf, open(output_file, 'w') as outfile:
            # Write header to output file
            outfile.write(f"# SNPs for {gene_name} ({chrom}:{start}-{end})\n")
            outfile.write("# CHROM\tPOS\tREF\tALT\tQUAL\tAA_CHANGE\n")
            
            for line in vcf:
                if line.startswith('#'):
                    continue
                
                # Parse VCF line
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                    
                chrom_from_vcf, pos_str, _, ref, alt = parts[0:5]
                qual = parts[5] if len(parts) > 5 else "."
                pos = int(pos_str)
                
                # Check if it's a SNP (both ref and alt are single bases)
                if len(ref) == 1 and all(len(a) == 1 for a in alt.split(',')):
                    # For negative strand, complement the alleles
                    working_ref = ref
                    working_alt = alt
                    if strand == "-":
                        working_ref = complement_base(ref)
                        working_alt = complement_base(alt)
                    
                    # Calculate amino acid change
                    aa_change = get_amino_acid_change(
                        working_sequence, 
                        pos, 
                        working_ref, 
                        working_alt, 
                        start if strand == "+" else end
                    )
                    
                    outfile.write(f"{chrom_from_vcf}\t{pos}\t{ref}\t{alt}\t{qual}\t{aa_change}\n")
                    snp_count += 1
        
        # Clean up
        if os.path.exists(temp_vcf):
            os.remove(temp_vcf)
            
        print(f"Extracted {snp_count} SNPs for {gene_name}")

def complement_base(base):
    """Return the complement of a nucleotide base."""
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return complements.get(base, base)

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return ''.join(complement_base(base) for base in reversed(sequence))

def main():
    # VCF file to process
    vcf_file = "output.vcf.gz"
    
    # Check if VCF file exists
    if not os.path.exists(vcf_file):
        print(f"Error: {vcf_file} not found.")
        sys.exit(1)
    
    # Ensure VCF file is indexed
    if not ensure_vcf_index(vcf_file):
        print("Failed to index VCF file. Using slower region extraction.")
    
    # Check chromosome naming in VCF
    check_cmd = "bcftools view -h output.vcf.gz | grep contig"
    contig_info = run_command(check_cmd)
    
    # Determine if chromosomes use "chr" prefix
    use_chr_prefix = False
    if contig_info and "chr" in contig_info:
        use_chr_prefix = True
    
    chr_prefix = "chr" if use_chr_prefix else ""
    
    # Try to determine if VCF is using hg19 or hg38
    # This is a simplistic approach - might need refinement
    check_variant_cmd = "bcftools view output.vcf.gz | head -n 100"
    sample_variants = run_command(check_variant_cmd)
    
    use_hg38 = True  # Default to hg38
    if sample_variants:
        # Simple heuristic - you might want to improve this
        if "7:34694" in sample_variants or "10:1158043" in sample_variants:
            use_hg38 = False
            print("Detected hg19/GRCh37 coordinates in VCF")
        else:
            print("Assuming hg38/GRCh38 coordinates in VCF")
    
    # Extract SNPs with amino acid changes
    extract_snps_with_aa_changes(vcf_file, use_hg38, chr_prefix)
    
    print("\nAll genes processed. SNPs extracted to gene_variants directory.")

if __name__ == "__main__":
    main()