#!/bin/bash -e
#SBATCH --job-name=herad_funannotate # job name (shows up in the queue)
#SBATCH --time=08:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=32 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
#SBATCH --dependency=afterok:49087996 # wait till previous purge haplotype task is complete

module purge
module load Apptainer
module load GeneMark-ES
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/11_funannotate_predict

export GENEMARK_PATH=/opt/nesi/CS400_centos7_bdw/GeneMark-ES/4.71-GCC-11.3.0/

#cp ../04_funannotate/output/DH641_Ppo_masked.fasta DH641_Ppo.fasta

#sed 's/.1_RagTag//' DH641_Ppo.fasta > DH641_Ppo_modified.fasta

# remove bad contig identified in a previous run
#awk 'BEGIN {p=1} /^>/{if ($0 ~ /^>k141_23843(\s|\r)?$/) p=0; else p=1} p' DH641_Ppo_modified.fasta > DH641_Ppo_badcontiqrmv.fasta

apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate predict -i ./00_purge_haplotigs/curated.fasta -o haploidDH641 \
    --species "Pleurotus purpureo-olivaceus" --strain DH641_haploid \
    --busco_seed_species coprinus --cpus 32


