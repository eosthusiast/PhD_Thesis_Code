#!/bin/bash -e
#SBATCH --job-name=herad_kraken2trial # job name (shows up in the queue)
#SBATCH --time=30:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=500GB          # Memory in MB
# #SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module purge
module load Kraken2
#module load aria2/1.36.0-GCCcore-11.3.0

cd /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L01_first_run_20240212/trimmed/Ppo/Assemblies_spades/
#cd /nesi/nobackup/ga03488/David_H/kraken2customDB/bigFungiDB

kraken2 --db /nesi/nobackup/ga03488/David_H/kraken2customDB/bigFungiDB/bigDB/ --output ICMP9630_kraken2_output_bigDBcust.txt --report ICMP9630_kraken2_report_bigDBcust.txt ICMP9630_Ppo_UDB-108_1_spades.fasta
kraken2 --db /nesi/nobackup/ga03488/David_H/kraken2customDB/bigFungiDB/bigDB/ --output DH639_kraken2_output_bigDBcust.txt --report DH639_kraken2_report_bigDBcust.txt DH639_Ppo_UDB-114_1_spades.fasta

#cd /nesi/project/ga03488/David/05_RPB2_checks/ref/Pleu

#for file in *.fna
#do
#    kraken2-build --add-to-library $file --db /nesi/nobackup/ga03488/David_H/kraken2customDB/bigFungiDB/bigDB/
#done

#kraken2-build --build --db /nesi/nobackup/ga03488/David_H/kraken2customDB/bigFungiDB/bigDB/ --threads 16