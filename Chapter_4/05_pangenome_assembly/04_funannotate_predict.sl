#!/bin/bash -e
#SBATCH --job-name=herad_funannotate # job name (shows up in the queue)
#SBATCH --time=22:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49104870 # wait till previous purge haplotype task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

#!/bin/bash -e
module purge
module load Apptainer
module load GeneMark-ES
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

export GENEMARK_PATH=/opt/nesi/CS400_centos7_bdw/GeneMark-ES/4.71-GCC-11.3.0/

# Loop over each .fasta file in the curated_genomes directory
for fasta_file in ../01_purged_noragtag/curated_genomes/*.fasta; do
    # Extract the basename (without extension)
    basename=$(basename "$fasta_file" .fasta)
    
    # Define the input and output files
    input_file="$fasta_file"
    output_file="../01_purged_noragtag/curated_genomes/${basename}_renamed.fasta"
    
    # Create the output file or clear it if it already exists
    > "$output_file"
    
    # Rename the fasta file headers after phasing
    while IFS= read -r line; do
        if [[ $line == \>* ]]; then
            # If the line starts with ">", keep only the part before the first space
            echo "${line%% *}" >> "$output_file"
        else
            # Otherwise, just write the line as is
            echo "$line" >> "$output_file"
        fi
    done < "$input_file"
    
    echo "Renaming completed. Output written to $output_file"
    
    # Remove bad contig identified in a previous run
    awk 'BEGIN {p=1} /^>/{if ($0 ~ /^>k141_23843(\s|\r)?$/) p=0; else p=1} p' "$output_file" > "${basename}_renamed.fasta"
    awk 'BEGIN {p=1} /^>/{if ($0 ~ /^>k141_19127(\s|\r)?$/) p=0; else p=1} p' "${basename}_renamed.fasta" > "${basename}_renamed2.fasta"
    mv "${basename}_renamed2.fasta" "${basename}_renamed.fasta"
    
    apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate mask -i "${basename}_renamed.fasta" -o "${basename}_masked.fasta"
    
    # Check if the directory exists, and create it if it does not
    if [ ! -d "${basename}_haploid" ]; then
        mkdir "${basename}_haploid"
    else
        echo "Directory ${basename}_haploid already exists. Skipping directory creation."
    fi
    
    # Run funannotate predict
    apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate predict -i "${basename}_masked.fasta" -o "${basename}_haploid" \
        --species "Pleurotus purpureo-olivaceus" --strain "${basename}_haploid" \
        -p /nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/haploidDH641_noragtag/predict_results/pleurotus_purpureo-olivaceus_dh641_haploid_noragtag.parameters.json --cpus 16
done



