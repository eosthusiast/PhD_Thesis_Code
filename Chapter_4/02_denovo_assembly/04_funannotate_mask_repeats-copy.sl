#!/bin/bash -e
#SBATCH --job-name=herad_funannotate # job name (shows up in the queue)
#SBATCH --time=01:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32768MB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module purge
module load Apptainer
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/04_funannotate/otherPleurotus/

#apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate mask -i DH639_Ppo_megahit_ragtag.fasta --cpus 16 -o DH639_masked.fasta

for file in *.fasta
do
    # Extract the basename (without extension) of the .fasta file
    basename=$(echo ${file} | cut -d'_' -f1-2)
    
    apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate mask -i $file --cpus 16 -o "./output/${basename}_masked.fasta"
done