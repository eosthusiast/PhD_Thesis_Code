#!/bin/bash -e
#SBATCH --job-name=funfinder # job name (shows up in the queue)
#SBATCH --time=2:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:50038244 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module load Python

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/05_FunFinder_Pangenome


python ../FunFinder_Pangenome.py -d ./input/ -o output -c 16 --EFFECTORP_PATH /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/EffectorP/EffectorP-3.0/EffectorP.py -a phobius, merops, dbcan, pfam, iprscan, antismash, tmhmm, signalp, effectors -p 0.5

# python ./EffectorP/EffectorP-3.0/EffectorP.py -f -o test.txt -i ./05_FunFinder_Pangenome/input/DH1310_haploid/DH1310_proteins.fasta
