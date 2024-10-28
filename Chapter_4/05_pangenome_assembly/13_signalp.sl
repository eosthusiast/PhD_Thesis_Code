#!/bin/bash -e
#SBATCH --job-name=signalp # job name (shows up in the queue)
#SBATCH --time=15:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=16GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:50038244 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
module load SignalP

#signalp6 --fastafile /path/to/input.fasta --organism eukarya --output_dir path/to/be/saved --format txt --mode fast

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

# Loop through all gbk files in the current directory
for fasta in ./*/predict_results/*proteins.fa; do
    #Extract the basename of the first folder (i.e., the folder containing the current file), excluding ./
    first_folder=$(basename $(dirname $(dirname "$fasta")))   
    
    echo "starting SignalP run on $fasta now"
    
    signalp6 --fastafile $fasta --organism eukarya --output_dir ./$first_folder/annotate_misc/signalp --format none --mode fast -wp 16

    echo "finished run on $fasta , moving on..."
    
done




