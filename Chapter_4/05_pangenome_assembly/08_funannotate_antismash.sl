#!/bin/bash -e
#SBATCH --job-name=antismash # job name (shows up in the queue)
#SBATCH --time=20:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49254984 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
module load antiSMASH
unset I_MPI_PMI_LIBRARY

#export PATH=/opt/nesi/CS400_centos7_bdw/antiSMASH/6.0.1-gimkl-2020a-Python-3.8.2/bin/meme:$PATH

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

# Loop through all gbk files in the current directory
for gbk in ./*/predict_results/*.gbk; do
    #Extract the basename of the first folder (i.e., the folder containing the current file), excluding ./
    first_folder=$(basename $(dirname $(dirname "$gbk")))
        
    outdir="./$first_folder/annotate_misc/antismash/"

    if [ ! -d "$outdir" ]; then
      mkdir -p "$outdir"
    fi      
    echo "starting annotation run on $first_folder now"
    antismash -c 16 --taxon fungi --output-dir $outdir --genefinding-tool glimmerhmm $gbk 
    echo "completed annotation run on $first_folder , moving on to next... "
done




