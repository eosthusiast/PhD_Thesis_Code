#!/bin/bash -e
#SBATCH --job-name=RepeatModeler # job name (shows up in the queue)
#SBATCH --time=2:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

# Load necessary modules
module purge
module load RepeatMasker

# Define paths
REF_GENOME=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/haploidDH641_noragtag/annotate_results/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa
OUTPUT_DIR=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/RepeatModeler2

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/RepeatModeler2

ln -s RM_*/consensi.fa.classified ./

RepeatMasker -pa 4 -e ncbi -gff -lib consensi.fa.classified -dir MaskerOutput $REF_GENOME