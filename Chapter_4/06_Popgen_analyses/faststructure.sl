#!/bin/bash -e
#SBATCH --job-name=faststructure # job name (shows up in the queue)
#SBATCH --time=2:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49254984 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/14_whole_genome_snpcalling/02_variants/plink_WG_SNPs
module purge
module load fastStructure

outdir=./structureselector

if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

for k in `seq 1 6`; do for i in `seq 1 10`; do structure.py -K $k --input=WG_SNPs --output=$outdir/out_run.$i --cv=10 \
--seed=$RANDOM; done; done
zip -j -q faststr_results.zip ./structureselector/*.meanQ *.log 

echo "Done, load into structureselector now."