#!/bin/bash -e
#SBATCH --job-name=merqury # job name (shows up in the queue)
#SBATCH --time=2:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module purge
module load Merqury
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/DH641_merqury

export MERQURY=/opt/nesi/CS400_centos7_bdw/Merqury/1.3-Miniconda3/merqury/
ln -s /opt/nesi/CS400_centos7_bdw/Merqury/1.3-Miniconda3/merqury/merqury.sh


# 1. Build meryl dbs on each input, could be run on different nodes
meryl k=17 count output read1.meryl /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L02_first_run_20240625/trimmed/DH641_combined_1.sub_1.fastq threads=16 memory=32g
meryl k=17 count output read2.meryl /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L02_first_run_20240625/trimmed/DH641_combined_1.sub_2.fastq threads=16 memory=32g

# 2. Merge
meryl union-sum output DH641.meryl read*.meryl threads=16 memory=32g

./merqury.sh DH641.meryl /nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/DH641_masked_noragtag.fasta merqury
