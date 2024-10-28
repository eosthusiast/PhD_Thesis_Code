#!/bin/bash -e
#SBATCH --job-name=pankmer_retrim # job name (shows up in the queue)
#SBATCH --time=0:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=150GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module load Python/3.11.6-foss-2023a
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/04_pankmer-ppofastq

#pankmer index -p /nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/01_trimmed/paired -o outPpoFastqretrim/ -t 16 --time
pankmer adj-matrix -i outPpoFastqretrim/ -o PpoFastqRetrim.csv
pankmer clustermap -i PpoFastqRetrim-24.csv -o PpoFastqRetrim_heatmap.svg --width 14 --height 14 --dend-ratio 0.35 --colormap viridis
pankmer collect -i outPpoFastqretrim/ -o PpoFastqRetrim_heatmap_collect.svg
pankmer count -i outPpoFastqretrim/

