#!/bin/bash

#SBATCH --job-name=Pu_tntblast_1       # Job name 
#SBATCH --time=01:23:13                # Walltime 
#SBATCH --mem=32768MB                  # Memory in MB
#SBATCH --account=ga03488              # Project code
#SBATCH --cpus-per-task=16             # CPUs per task 

# Work directory (adjust as needed)
#work_dir=/nesi/nobackup/ga03488/David_H/PPoutput/L01_first_run_20240212/trimmed/Pul
#work_dir=/nesi/nobackup/ga03488/David_H/PPoutput/L01_first_run_20240212/trimmed/Dja+Pars+Ppo_loci
work_dir=/nesi/project/ga03488/David/

cd $work_dir 

# Define folders
input_folder="Reference_genomes"
target_folder="work_dir=/nesi/nobackup/ga03488/David_H/PPoutput/L01_first_run_20240212/trimmed/Pul/Target_gene_regions"
output_folder="locus_extraction_0503"

# Create output folder if it doesn't exist
mkdir -p $output_folder

# Iterate through input and target files
for d_file in $input_folder/*.fasta; do
    d_basename=$(basename $d_file "-d")  # Extract base filename
    d_prefix=$(echo $d_basename | cut -d'_' -f1-2)  # Get first two parts before '_'

    for i_file in $target_folder/*.txt; do
        i_target=$(basename $i_file "-i" | cut -d'_' -f1)  # Get target gene
        
                # Print paths for inspection 
        echo "Input file: $d_file"
        echo "Target file: $i_file"

        output_file="${output_folder}/${d_prefix}_${i_target}_tntblast.fasta"

        # Run tntblast and extract first sequence
        /nesi/project/ga03488/software/tntblast/thermonucleotideBLAST-2.61/tntblast -e 40 -E 45 -i $i_file -d $d_file -o ${d_prefix}_${i_target}_tmp.fasta -m 1 -l 3000 
        head -2 ${d_prefix}_${i_target}_tmp.fasta | sed "s/>.*/>${d_prefix}_${i_target}/" > $output_file

    done
done