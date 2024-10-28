#!/bin/bash -e
#SBATCH --job-name=RepeatModeler # job name (shows up in the queue)
#SBATCH --time=6:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

# Load necessary modules
module purge
module load RepeatModeler
module load trf/4.09

# Define paths
REF_GENOME=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/haploidDH641_noragtag/annotate_results/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa
OUTPUT_DIR=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/RepeatModeler2

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/RepeatModeler2

# Build database
BuildDatabase -name DH641_Ppo -engine ncbi ../haploidDH641_noragtag/annotate_results/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa

# Run RepeatModeler2
RepeatModeler -database DH641_Ppo -engine ncbi -pa 16 -LTRStruct > out.log