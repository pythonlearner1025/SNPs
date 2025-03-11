#!/bin/bash

# Script to download gnomAD exome VCF files for chromosomes 1-22
# The files will be downloaded to the current working directory
# If script is run again, it resumes from the highest chromosome number
# that exists in the directory (including partial downloads)

# Base URL for gnomAD exome VCF files
BASE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes"

# Create a log file
LOG_FILE="download_log.txt"
echo "Starting download at $(date)" | tee -a $LOG_FILE

# Find the highest chromosome number that exists (complete or partial)
HIGHEST_CHR=0
for chr in {1..22}; do
    FILENAME="gnomad.exomes.v4.1.sites.chr$chr.vcf.bgz"
    if [ -f "$FILENAME" ]; then
        HIGHEST_CHR=$chr
    fi
done

# If we found existing files, start from the highest one
# Otherwise start from chromosome 1
START_CHR=1
if [ $HIGHEST_CHR -gt 0 ]; then
    START_CHR=$HIGHEST_CHR
    echo "Resuming downloads from chromosome $START_CHR (inclusive)" | tee -a $LOG_FILE
else
    echo "Starting fresh download from chromosome 1" | tee -a $LOG_FILE
fi

# Download VCF files from START_CHR to chromosome 22
for chr in $(seq $START_CHR 22); do
    echo "Downloading chromosome $chr..."
    
    # Construct the file URL
    FILE_URL="$BASE_URL/gnomad.exomes.v4.1.sites.chr$chr.vcf.bgz"
    
    # Construct the local filename
    FILENAME="gnomad.exomes.v4.1.sites.chr$chr.vcf.bgz"
    
    # Download the file using wget
    # -c enables resume of interrupted downloads
    wget -c $FILE_URL -O $FILENAME
    
    # Check if download was successful
    if [ $? -eq 0 ]; then
        echo "Successfully downloaded $FILENAME" | tee -a $LOG_FILE
    else
        echo "Failed to download $FILENAME" | tee -a $LOG_FILE
    fi
done

echo "All downloads completed at $(date)" | tee -a $LOG_FILE
echo "Files have been saved to $(pwd)"

# Optional: Verify the downloaded files
echo "Verifying downloaded files..."
for chr in {1..22}; do
    FILENAME="gnomad.exomes.v4.1.sites.chr$chr.vcf.bgz"
    if [ -f "$FILENAME" ]; then
        SIZE=$(du -h "$FILENAME" | cut -f1)
        echo "$FILENAME exists ($SIZE)" | tee -a $LOG_FILE
    else
        echo "$FILENAME is missing!" | tee -a $LOG_FILE
    fi
done
