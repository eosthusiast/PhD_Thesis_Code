#!/bin/bash -e
#SBATCH --job-name=orthofinder # job name (shows up in the queue)
#SBATCH --time=2:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=72 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:50038244 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
module load OrthoFinder
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes


echo #Starting orthofinder run now"

orthofinder -f ../04_orthofinder/ -S diamond_ultra_sens -t 72

