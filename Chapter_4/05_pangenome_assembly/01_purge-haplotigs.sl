#!/bin/bash -e
#SBATCH --job-name=01_haplotigs1 # job name (shows up in the queue)
#SBATCH --time=18:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module purge
module load BWA
module load SAMtools
module load purge_haplotigs
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/01_purged_noragtag

# Directory paths
genome_dir="/nesi/nobackup/ga03488/David_H/Ppo_Popgen/03_ragtag_tuber-regium"
fastq_dir="/nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L01_first_run_20240212/trimmed/Ppo/trimmed_fastq"

# Iterate over each genome file
for genome_file in "$genome_dir"/*.fasta; do
    # Extract the basename (part before the first '_')
    basename=$(basename "$genome_file" | cut -d'_' -f1)

    # Construct the full paths for the fastq files
    fastq1_file=$(find "$fastq_dir" -name "${basename}*_1.trim_1.fastq.gz" -print -quit)
    fastq2_file=$(find "$fastq_dir" -name "${basename}*_1.trim_2.fastq.gz" -print -quit)

    # Replace DH640 with the current basename
    output_bam="${basename}_selfmap.bam"

    # Print paths to verify (optional, can be removed)
    echo "Genome: $genome_file"
    echo "Fastq 1: $fastq1_file"
    echo "Fastq 2: $fastq2_file"
    echo "Output BAM: $output_bam"

    # Run BWA index
    bwa index "$genome_file"

    # Map and sort
    bwa mem "$genome_file" -t 16 "$fastq1_file" "$fastq2_file" | samtools sort -@ 15 -m 1G -o "$output_bam" -T ali.tmp

    # Get stats on mapping
    samtools flagstat "$output_bam"

    # Generate a read-depth histogram
    purge_haplotigs readhist -b "$output_bam" -g "$genome_file" -t 16 -d 200
    
    #
    echo "$genome_file complete, moving on to next sample "

done

# Next step: inspect and find manual cutoffs from histograms, then configure and run 02_purge-haplotigs.sl
