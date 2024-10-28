#!/bin/bash -e
#SBATCH --job-name=herad_QUAST_1 # job name (shows up in the queue)
#SBATCH --time=8:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32768MB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

# cd /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L01_first_run_20240212/trimmed/Ppo/megahit/Assemblies
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/01_Draft_assemblies/QUAST

module load QUAST

for FILENAME in *.fasta; do
    # Extract first string before "_"
    folder_name="${FILENAME%%_*}"

    # Create output folder if it doesn't exist
    mkdir -p "./${folder_name}"

    # Run QUAST command with the generated folder
    quast.py "$FILENAME" --split-scaffolds -t 16 --fungus -o "../QUAST/${folder_name}"
done