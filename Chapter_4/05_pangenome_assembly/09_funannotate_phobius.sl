#!/bin/bash -e
#SBATCH --job-name=phobius # job name (shows up in the queue)
#SBATCH --time=0:30:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49254984 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
#module load antiSMASH

#export PATH=/opt/nesi/CS400_centos7_bdw/antiSMASH/6.0.1-gimkl-2020a-Python-3.8.2/bin/meme:$PATH

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

# Loop through all gbk files in the current directory
#for fasta in ./*/predict_results/*proteins.fa; do
#    #Extract the basename of the first folder (i.e., the folder containing the current file), excluding ./
#    first_folder=$(basename $(dirname $(dirname "$fasta")))   
#    
#    echo "starting phobius run on $fasta now"
#    perl /nesi/project/ga03488/software/phobius/phobius.pl $fasta -short > ./$first_folder/annotate_misc/phobius.results.txt
#    echo "finished phobius run on $fasta"
#    
#done


perl /nesi/project/ga03488/software/phobius/phobius.pl /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes/DH1310_haploid/predict_results/Pleurotus_purpureo-olivaceus_DH1310_haploid.proteins.fa -short > /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes/DH1310_haploid/annotate_misc/phobius_result.txt