#!/bin/bash -e
#SBATCH --job-name=herad_kraken2decontam # job name (shows up in the queue)
#SBATCH --time=8:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

# Load any necessary modules (if required)
# module load python

# Navigate to the directory containing the .fasta files
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/01_Draft_assemblies/

# Loop through each .fasta file in the current directory
for fasta_file in *.fasta
do
    # Extract the basename (without extension) of the .fasta file
    base_name=$(basename "$fasta_file" .fasta)
    
    # Construct the output filenames based on the basename
    out_file="${base_name}_out.txt"
    report_file="${base_name}_report.txt"
    decontam_file="${base_name}_decontam.fasta"
    
    # Run the Python script with the appropriate arguments
    python /nesi/nobackup/ga03488/David_H/kraken2customDB/extract_kraken_reads.py \
        -k "$out_file" \
        -r "$report_file" \
        -s "$fasta_file" \
        -o "$decontam_file" \
        --taxid 2 4890 554915 33208 33090 10239 \
        --exclude \
        --include-children
done
