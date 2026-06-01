#!/usr/bin/bash

# Set base URLs
FASTA_URL="https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gtf.gz"

mkdir -p references

# Download the FASTA chr1 file
curl -L "${FASTA_URL}" -o references/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

# Download the GTF file
curl -L "${GTF_URL}" -o references/Homo_sapiens.GRCh38.114.chr.gtf.gz

# Unzip the files

gunzip references/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
gunzip references/Homo_sapiens.GRCh38.114.chr.gtf.gz

# separate only chr1 annotations from the GTF file
awk '$1 == "1"' references/Homo_sapiens.GRCh38.114.chr.gtf > references/Homo_sapiens.GRCh38.114.chr1.gtf

# remove complete gtf file
rm references/Homo_sapiens.GRCh38.114.chr.gtf

echo -e "\nHuman reference chromosome 1 - genome and annotation files downloaded.\n"

# Download the fastq files

mkdir -p samples

FASTQ_URL="https://www.encodeproject.org/files/"

curl -L "${FASTQ_URL}"ENCFF309DAU/@@download/ENCFF309DAU.fastq.gz -o samples/ENCFF309DAU.fastq.gz
curl -L "${FASTQ_URL}"ENCFF168OKB/@@download/ENCFF168OKB.fastq.gz -o samples/ENCFF168OKB.fastq.gz


echo -e "\nENCODE Project CapTrap cDNA - Downloaded ENCFF309DAU and ENCFF168OKB fastq files for heart samples.\n"