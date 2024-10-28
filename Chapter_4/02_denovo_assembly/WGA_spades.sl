#!/bin/bash -e
#SBATCH --job-name=herad_spades_1 # job name (shows up in the queue)
#SBATCH --time=40:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32768MB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

cd /nesi/nobackup/ga03488/David_H/PPoutput/L01_first_run_20240212/trimmed/Ppo

bash spadesWGA_script.sh