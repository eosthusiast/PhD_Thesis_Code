#!/bin/bash
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/03_variant_calling/

module load BCFtools

VCF_IN=raw_variants_renamed.vcf
VCF_OUT=vcffiltered_SNPs.vcf.gz

# Filter for SNPs and only bi-allelic sites
bcftools view --threads 4 -m 2 -M 2 -v snps \
    -O z -o biallelic.vcf.gz $VCF_IN

# Annotate SNPset with HWE and MAF annotations
bcftools +fill-tags biallelic.vcf.gz -O z -o biallelic_annotated.vcf.gz -- -t HWE,MAF

VCF_IN=masked_variants.vcf
OUT=./total/global_variant

module load VCFtools

# from https://speciationgenomics.github.io/variant_calling/
#Calculate allele frequency
vcftools --vcf $VCF_IN --freq2 --out $OUT
#Calculate mean depth per individual
vcftools --vcf $VCF_IN --depth --out $OUT
#Calculate mean depth per site
vcftools --vcf $VCF_IN --site-mean-depth --out $OUT
#Calculate site quality
vcftools --vcf $VCF_IN --site-quality --out $OUT
#Calculate proportion of missing data per individual
vcftools --vcf $VCF_IN --missing-indv --out $OUT
#Calculate proportion of missing data per site
vcftools --vcf $VCF_IN --missing-site --out $OUT
#Calculate heterozygosity and inbreeding coefficient per individual
vcftools --vcf $VCF_IN --het --out $OUT

# run R script VariantGraphing.R

# set filters
MAF=0.1
MISS=1
QUAL=30
MIN_DEPTH=20
MAX_DEPTH=150

# perform the filtering with vcftools - no MAF and no Missingness
vcftools --gzvcf $VCF_IN \
--minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

# perform the filtering with vcftools
#vcftools --gzvcf $VCF_IN \
#--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
#--minQ $QUAL \
#--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
#--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
#$VCF_OUT