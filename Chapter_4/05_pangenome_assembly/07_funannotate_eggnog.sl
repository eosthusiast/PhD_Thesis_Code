#!/bin/bash -e
#SBATCH --job-name=eggnog # job name (shows up in the queue)
#SBATCH --time=24:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49254984 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
module load eggnog-mapper/2.1.12-gimkl-2022a

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

export EGGNOG_DATA_DIR=/nesi/nobackup/ga03488/David_H/funannotate/eggnogdb/

# Loop through all .proteins.fa files in the current directory
for proteins in ./*/predict_results/*.proteins.fa; do
    #Extract the basename of the first folder (i.e., the folder containing the current file), excluding ./
    first_folder=$(basename $(dirname $(dirname "$proteins")))
        
    # Run the emapper.py command
    emapper.py -i "$proteins" -o "./$first_folder/$first_folder" --tax_scope 4751 --cpu 16 --override
        
    mv ./$first_folder/*.annotations ./$first_folder/annotate_misc/eggnog.emapper.annotations
    
done




