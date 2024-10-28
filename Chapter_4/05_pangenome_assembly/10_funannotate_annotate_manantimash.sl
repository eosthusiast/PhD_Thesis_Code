#!/bin/bash -e
#SBATCH --job-name=funannotate # job name (shows up in the queue)
#SBATCH --time=2:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:50028563 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
module load Apptainer
module load SignalP
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

# Loop over each .fasta file in the curated_genomes directory
for fasta in ./*/predict_results/*proteins.fa; do
    #Extract the basename of the first folder (i.e., the folder containing the current file), excluding ./
    first_folder=$(basename $(dirname $(dirname "$fasta")))   
echo "starting annotation run on $fasta now"
    apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate annotate -i $first_folder --cpus 16 --phobius ./"$first_folder"/annotate_misc/phobius_result.txt --antismash ./"$first_folder"/annotate_misc/antismash/*_haploid.gbk
echo "completed annotation run on $fasta , moving on..."
done


