#!/bin/bash -e
#SBATCH --job-name=01_retrim # job name (shows up in the queue)
#SBATCH --time=6:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module load cutadapt
module load TrimGalore

# First, adapter and end trimming using MGI adapters, without /5Phos/ start of adapter1

raw=/nesi/nobackup/ga03488/David_H/refreshed_backup/PE150_deconvoluted/L01_renamed/Ppo/
trim_out=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/01_trimmed
reports=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/01_trimmed/reports

for samp in ${raw}*_1.fq.gz
do
base=$(basename ${samp} _1.fq.gz)
echo "Running Trim_galore for ${base}..."
trim_galore --paired \
    --adapter AGTCGGAGGCCAAGCGGTCTTAGGAAGACAATCAG \
    --adapter2 TTGTCTTCCTAAGCAACTCCTTGGCTCACAGAACGACATGGCTACGATCCGACTT \
    --quality 30 \
    --cores 8 \
    --fastqc_args "--nogroup --outdir ${reports} --threads 16" \
    --length 50 \
    --output_dir ${trim_out}/ \
    --clip_R1 11 \
    --clip_R2 11 \
    --three_prime_clip_R1 4 \
    --three_prime_clip_R2 4 \
    --retain_unpaired \
    --length_1 55 \
    --length_2 55 \
    ${raw}${base}_1.fq.gz \
    ${raw}${base}_2.fq.gz
mv ${trim_out}/*.txt ${reports}
done
