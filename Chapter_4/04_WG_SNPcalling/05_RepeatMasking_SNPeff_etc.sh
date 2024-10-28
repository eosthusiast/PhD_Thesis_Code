#!/bin/bash
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/03_variant_calling/
VCF_IN=masked_variants.vcf


## Now continuing with removing repeats from SNPset
# Convert RepeatMasker output to BED file for mapping to VCF snpset
# https://github.com/4ureliek/Parsing-RepeatMasker-Outputs/blob/master/RMout_to_bed.pl
ref=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/RepeatModeler2/MaskerOutput/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa.out
perl RMout_to_bed.pl $ref
bed=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/RepeatModeler2/MaskerOutput/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa.out.bed
module load BEDTools
# Extract the header from the original VCF
tabix -H biallelic.vcf.gz > header.txt
bgzip -d -k biallelic.vcf.gz
bedtools intersect -v -a biallelic.vcf  -b $bed > filtered_variants.vcf
cat header.txt filtered_variants.vcf > masked_variants.vcf
rm header.txt filtered_variants.vcf
# Count variants in the filtered VCF file
bcftools view -H masked_variants.vcf | wc -l # 822534 SNPs

# After that I did more filtering with bcftools to remove sites with missing genotype information
bcftools view --threads 4 -i 'N_PASS(GT=="mis")<1' \
  -O z -o bcffiltered_SNPs.vcf.gz \
  $VCF_IN
bcftools view -H bcffiltered_SNPs.vcf.gz | wc -l # 458886 SNPs

### SNPeff of all
module load snpEff
java -Xmx15g -jar $EBROOTSNPEFF/snpEff.jar -c /nesi/project/ga03488/my_snpEff.config -v -ud 1000 -stats stats_1k_WG_all.html DH641 bcffiltered_SNPs.vcf.gz > variants_snpeff_1k_WG_all.vcf

# subset vcfs to filter introns, exons, etc
module load BCFtools
bcftools view -i 'INFO/ANN[*] ~ "missense_variant" || INFO/ANN[*] ~ "synonymous_variant" || INFO/ANN[*] ~ "stop_gained" || INFO/ANN[*] ~ "stop_lost" || INFO/ANN[*] ~ "start_lost"' variants_snpeff_1k_WG_all.vcf -o coding_regions_only.vcf
bcftools view -i 'INFO/ANN[*] ~ "intron_variant" || INFO/ANN[*] ~ "UTR_variant" || INFO/ANN[*] ~ "non_coding_transcript_variant" || INFO/ANN[*] ~ "splice_region_variant" || INFO/ANN[*] ~ "regulatory_region_variant"' variants_snpeff_1k_WG_all.vcf -o noncoding_regions_only.vcf
bcftools view -i 'INFO/ANN[*] ~ "upstream_gene_variant" || INFO/ANN[*] ~ "downstream_gene_variant"' variants_snpeff_1k_WG_all.vcf -o up_downstream_regions.vcf
bcftools view -i 'INFO/ANN[*] ~ "intergenic_region"' variants_snpeff_1k_WG_all.vcf -o intergenic_regions.vcf

#80% subset
bcftools view -i 'COUNT(GT="het" | GT="alt") >= 22' bcffiltered_SNPs.vcf.gz -o shared_snps_80p.vcf
bcftools view -H shared_snps_80p.vcf | wc -l ## 62,700 SNPs 

#96% subset
bcftools view -i 'COUNT(GT="het" | GT="alt") >= 27' bcffiltered_SNPs.vcf.gz -o shared_snps_96p.vcf
bcftools view -H shared_snps_96p.vcf | wc -l ## 61,904 SNPs 

#100% subset
bcftools view -i 'COUNT(GT="het" | GT="alt") >= 28' bcffiltered_SNPs.vcf.gz -o shared_snps_100p.vcf
bcftools view -H shared_snps_100p.vcf | wc -l ## 61,904 SNPs 

### SNPeff of 80%
module load snpEff
java -Xmx15g -jar $EBROOTSNPEFF/snpEff.jar -c /nesi/project/ga03488/my_snpEff.config -v -ud 1000 -stats stats_1k_WG_80.html DH641 shared_snps_80p.vcf > variants_snpeff_1k_WG_80.vcf

# subset vcfs to filter introns, exons, etc
module load BCFtools
bcftools view -i 'INFO/ANN[*] ~ "missense_variant" || INFO/ANN[*] ~ "synonymous_variant" || INFO/ANN[*] ~ "stop_gained" || INFO/ANN[*] ~ "stop_lost" || INFO/ANN[*] ~ "start_lost"' variants_snpeff_1k_WG_80.vcf -o coding_regions_only.vcf
bcftools view -i 'INFO/ANN[*] ~ "intron_variant" || INFO/ANN[*] ~ "UTR_variant" || INFO/ANN[*] ~ "non_coding_transcript_variant" || INFO/ANN[*] ~ "splice_region_variant" || INFO/ANN[*] ~ "regulatory_region_variant"' variants_snpeff_1k_WG_80.vcf -o noncoding_regions_only.vcf
bcftools view -i 'INFO/ANN[*] ~ "upstream_gene_variant" || INFO/ANN[*] ~ "downstream_gene_variant"' variants_snpeff_1k_WG_80.vcf -o up_downstream_regions.vcf
bcftools view -i 'INFO/ANN[*] ~ "intergenic_region"' variants_snpeff_1k_WG_80.vcf -o intergenic_regions.vcf