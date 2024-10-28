#!/bin/bash -e
#SBATCH --job-name=herad_haplotigs1 # job name (shows up in the queue)
#SBATCH --time=4:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=32 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module purge
module load BWA
module load SAMtools
module load purge_haplotigs
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/11_funannotate_predict/00_purge_haplotigs

# Map and sort
#bwa mem ../trialDH641/annotate_results/Pleurotus_purpureo-olivaceus_DH641.scaffolds.fa -t 32 \
# /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L02_first_run_20240625/trimmed/DH641_combined_1.trim_1.fastq.gz #/nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L02_first_run_20240625/trimmed/DH641_combined_1.trim_2.fastq.gz \
# | samtools sort -@ 8 -m 1G -o PH.bam -T ali.tmp

# Get stats on mapping
#samtools flagstat PH.bam

# Generate a read-depth histogram to determine SNP filtering cutoffs
#purge_haplotigs readhist -b PH.bam -g ../trialDH641/annotate_results/Pleurotus_purpureo-olivaceus_DH641.scaffolds.fa -t 32 -d 7000

# Filter based on manual cutoffs
purge_haplotigs  cov  -i PH.bam.gencov  -l 150  -m 650  -h 1725

# Identify and remove haplotigs
purge_haplotigs purge  -g ../trialDH641/annotate_results/Pleurotus_purpureo-olivaceus_DH641.scaffolds.fa  -c coverage_stats.csv 

#purge_haplotigs purge  -g ../trialDH641/annotate_results/Pleurotus_purpureo-olivaceus_DH641.scaffolds.fa  -c coverage_stats.csv  -a 60 -o curated_a60

# Generate dotplots from final alignments - option: -a / -align_cov     Percent cutoff for identifying a contig as a haplotig. DEFAULT = 70
purge_haplotigs purge  -g ../trialDH641/annotate_results/Pleurotus_purpureo-olivaceus_DH641.scaffolds.fa  -c coverage_stats.csv  -a 70  -d  -b PH.bam

