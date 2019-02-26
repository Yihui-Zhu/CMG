# CMG

This script is used to generate SNPs located in Differential Methylated Regions or windows filesand
test for in cis regulation associations.

The R package is still under developing now. The R code is available.

## Setting up input files: BED files and VCF files
Differential Methylated Regions or windows are stored in bed file.
BED file includes chrom, chromStart, chromEnd
The format of BED format can be found at: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

VCF files are generated from WGS data

## Generate SNP count files
Use VCFtools for each vcf files to generate SNP in each window counts for each sample
To download and use VCFtools, please follow: https://vcftools.github.io/index.html

For each sample: 
vcftools --gzvcf input.vcf --bed input.bed --counts --out output_snp_count.csv

The output files are ending with frq.count
Put all samples frq.count files and bed file into the same directory

## Three Major Functions
Function1: Make summarize genotype data files for all samples.
Function2: Select the SNP located inside DMRs or windows.
Function3: Number the major homozygous genotype: 1; hetezygous genotype: 2; minor homozygous genotype: 3.
