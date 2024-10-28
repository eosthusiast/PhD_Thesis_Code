#!/bin/bash

#SBATCH --job-name=blastn_subject_1       # Job name 
#SBATCH --time=01:23:13                # Walltime 
#SBATCH --mem=32768MB                  # Memory in MB
#SBATCH --account=ga03488              # Project code
#SBATCH --cpus-per-task=16             # CPUs per task 

module load BLAST

# Work directory (adjust as needed)
#work_dir=/nesi/nobackup/ga03488/David_H/PPoutput/L01_first_run_20240212/trimmed/Pul
work_dir=/nesi/nobackup/ga03488/David_H/PPoutput/L01_first_run_20240212/trimmed/Dja+Pars+Ppo_loci

cd $work_dir 

# Define folders
input_folder="Assemblies_spades"
target_folder="Target_gene_regions/reference_seq"
output_folder="locus_extraction_blast"

# Create output folder if it doesn't exist
mkdir -p $output_folder

# Iterate through input and target files
for d_file in $input_folder/*.fasta; do
    d_basename=$(basename $d_file "-d")  # Extract base filename
    d_prefix=$(echo $d_basename | cut -d'_' -f1-2)  # Get first two parts before '_'

    for i_file in $target_folder/*.fasta; do
        i_target=$(basename $i_file "-i" | cut -d'_' -f1)  # Get target gene
        
                # Print paths for inspection 
        echo "Input file: $d_file"
        echo "Target file: $i_file"

        output_file="${output_folder}/${d_prefix}_${i_target}_blast_hit.fasta"

        # Run blast and extract first sequence
        
        blastn -query $i_file -subject $d_file -max_target_seqs 1 -outfmt '6 sseq' | head -n 1 | sed -e "s/^/>${d_prefix}_${i_target}_blast_hit\n/" > $output_file


    done
done