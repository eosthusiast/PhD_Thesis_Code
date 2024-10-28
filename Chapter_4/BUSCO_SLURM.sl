#!/bin/bash -e
#SBATCH --job-name=herad_BUSCO_1 # job name (shows up in the queue)
#SBATCH --time=8:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32768MB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module load BUSCO

cd /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L01_first_run_20240212/trimmed/Ppo/megahit/Assemblies

export NUMEXPR_MAX_THREADS=16

# Iterate through each assembly
for assembly_file in *.fasta; do
    assembly_name=$(basename "$assembly_file" .fasta)  # Extract assembly name

    # Run BUSCO 
    busco -i "$assembly_file" -l agaricales_odb10 -m genome -c 16 -o "QC/BUSCO_$assembly_name"
done
