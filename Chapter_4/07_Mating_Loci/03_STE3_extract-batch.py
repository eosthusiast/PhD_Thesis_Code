from Bio import SeqIO
import pandas as pd
import os
import glob

# Define the directory containing protein FASTA files
SUBJECT_DIR = "/nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/05_FunFinder_Pangenome/input/*/*_proteins.fasta"

# Create the output directory if it doesn't exist
output_dir = "./matB_files_diploid"
os.makedirs(output_dir, exist_ok=True)
output_dir_protein =  "./matB_proteins"
os.makedirs(output_dir_protein, exist_ok=True)

# File to store all ProteinName hits
protein_hits_file = os.path.join(output_dir, "ProteinName_hits.txt")
with open(protein_hits_file, "w") as hits_file:
    hits_file.write("ProteinName\tChromosome\tStart\tEnd\n")  # Add headers

# File to store samples where STE3_1 was not found
ste3_not_found_file = os.path.join(output_dir, "STE3_1_not_found.txt")
with open(ste3_not_found_file, "w") as not_found_file:
    not_found_file.write("STE3_1_not_found\n")  # Add header

# Function to generate sample_short by removing specific parts from sample name
def get_sample_short(sample):
    sample_short = sample.replace("Pleurotus_purpureo-olivaceus_", "")
    sample_short = sample_short.split('_')[0]  # Remove anything after the first underscore after the prefix
    return sample_short

# Loop through all gff3 files in the directory
for gff_file in os.listdir("./genomes_annotated"):
    if gff_file.endswith(".gff3"):
        sample = gff_file.split(".")[0]
        sample_short = get_sample_short(sample)  # Generate the sample_short
        fasta_file = f"./genomes_annotated/{sample}.scaffolds.fa"

        print(f"Processing {sample}, Sample Short: {sample_short}")

        # Check if the FASTA file exists
        if not os.path.exists(fasta_file):
            print(f"Missing FASTA file for {sample}, skipping...")
            continue

        # Load the GFF file
        gff_path = f"./genomes_annotated/{gff_file}"
        gff_df = pd.read_csv(gff_path, sep="\t", header=None, comment='#')
        gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

        ste3_found = False  # Flag to check if STE3_1 is found

        # Find STE3_* proteins
        for index, row in gff_df.iterrows():
            if 'STE3_' in row['attributes']:
                protein_id = row['attributes'].split("ID=")[1].split(";")[0]
                # Extract STE3_ID from the "Name=" field in the attributes
                STE3_ID = row['attributes'].split("Name=")[1].split(";")[0]
                chromosome = row['seqid']
                start_pos = max(1, row['start'] - 5000)  # Ensure start doesn't go below 1
                end_pos = row['end'] + 5000

                # Use the extracted STE3_ID in the ProteinName
                ProteinName = f"{sample_short}_{protein_id}-T1_{STE3_ID}"
                print(f"Found STE3 protein: {ProteinName}, Chromosome: {chromosome}, Start: {start_pos}, End: {end_pos}")

                # Check if STE3_1 was found
                if 'STE3_1' in STE3_ID:
                    ste3_found = True

                # Save the ProteinName to the hits file
                with open(protein_hits_file, "a") as hits_file:
                    hits_file.write(f"{ProteinName}\t{chromosome}\t{start_pos}\t{end_pos}\n")

                # Save the extracted GFF annotations
                gff_filtered = gff_df[(gff_df['seqid'] == chromosome) &
                                      (gff_df['start'] >= start_pos) &
                                      (gff_df['end'] <= end_pos)]
                extracted_gff = os.path.join(output_dir, f"{ProteinName}_extracted.gff")
                gff_filtered.to_csv(extracted_gff, sep="\t", header=False, index=False)
                print(f"Extracted GFF annotations written to {extracted_gff}")

                # Read the FASTA file and extract the relevant region
                output_fasta = os.path.join(output_dir, f"{ProteinName}_extracted.fasta")
                with open(fasta_file, "r") as fasta_handle:
                    for record in SeqIO.parse(fasta_handle, "fasta"):
                        if record.id == chromosome:
                            # Extract the specific region and write it to a new FASTA file
                            sub_record = record[start_pos-1:end_pos]  # Convert to 0-based indexing
                            SeqIO.write(sub_record, output_fasta, "fasta")
                            print(f"Extracted sequence written to {output_fasta}")

        # If STE3_1 was not found, write to the not_found_file
        if not ste3_found:
            with open(ste3_not_found_file, "a") as not_found_file:
                not_found_file.write(f"{sample_short}\n")
            print(f"STE3_1 not found in {sample_short}, logged.")