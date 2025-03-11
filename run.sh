#!/bin/bash
# Run Short Sleep Analysis Script
# This script runs the short sleep analysis on the family VCF files

# Install required packages if needed
echo "Installing required packages..."
pip install pyvcf pandas tqdm pybedtools

# Run the analysis
echo "Running short sleep analysis..."
python snp.py \
  --mom marie.vcf.gz \
  --grandma lavinia.vcf.gz \
  --aunt tess.vcf.gz \
  --son grant.vcf.gz \
  --output short_sleep_candidates.tsv \
  --min-frequency 1e-5

echo "Analysis complete! Results saved to short_sleep_candidates.tsv"
