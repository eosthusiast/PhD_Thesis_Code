#!/bin/bash -e
#SBATCH --job-name=03_varcall # job name (shows up in the queue)
#SBATCH --time=3:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# # SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
# #SBATCH --mail-type=ALL

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/03_variant_calling


# Create PLINK files
module load PLINK/1.09b6.16
VCF=ArthursPass_snps_80p.vcf
# perform linkage pruning - i.e. identify prune sites
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --indep-pairwise 50 10 0.1 --out ArthursPass_snps_80p

# prune and create pca
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --extract ArthursPass_snps_80p.prune.in --make-bed --pca --out ArthursPass_snps_80p
    
# create distance matrix
plink --bfile ArthursPass_snps_80p --recode --allow-extra-chr
plink --file plink --genome --allow-extra-chr
plink --file plink --cluster --matrix --allow-extra-chr
plink --file plink --read-genome plink.genome --cluster --mds-plot 4 --allow-extra-chr