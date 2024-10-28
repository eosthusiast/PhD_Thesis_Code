#!/bin/bash -e
#SBATCH --job-name=herad_pankmer_1 # job name (shows up in the queue)
#SBATCH --time=4:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module load Python/3.11.6-foss-2023a
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/02_pankmer

pankmer index -g ../01_purged_noragtag/curated_genomes/input_genomes.tar -o out/ -t 16 --time
pankmer adj-matrix -i out/ -o PpoHapNoRag.csv
pankmer clustermap -i PpoHapNoRag.csv -o PpoHapNoRag_heatmap.svg --width 14 --height 14 --dend-ratio 0.35 --colormap viridis