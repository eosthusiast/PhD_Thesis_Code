from Bio import SeqIO
import pandas as pd
import os

# Create the output directory if it doesn't exist
output_dir = "./matA_files_diploid"
os.makedirs(output_dir, exist_ok=True)

# Load the matA_gff_table.txt file
gff_table = pd.read_csv("matA_gff_table.txt", sep="\t")

# Function to generate sample_short by removing specific parts from sample name
def get_sample_short(sample):
    sample_short = sample.replace("Pleurotus_purpureo-olivaceus_", "")
    sample_short = sample_short.split('_')[0]  # Remove anything after the first underscore after the prefix
    return sample_short

# Loop through each row (each sample) in the table
for index, row in gff_table.iterrows():
    sample = row['Sample']
    chromosome = row['Chromosome']
    start_pos = int(row['Start'])
    end_pos = int(row['End'])
    sample_short = get_sample_short(sample)  # Generate the sample_short

    print(f"Processing {sample}, Chromosome: {chromosome}, Start: {start_pos}, End: {end_pos}, Sample Short: {sample_short}")

    # Define file paths for the GFF and FASTA files based on the sample name
    gff_file = f"./genomes_annotated/{sample}.gff3"
    fasta_file = f"./genomes_annotated/{sample}.scaffolds.fa"

    # Check if the GFF and FASTA files exist
    if not os.path.exists(gff_file) or not os.path.exists(fasta_file):
        print(f"Missing files for {sample}, skipping...")
        continue

    # Load the GFF file
    gff_df = pd.read_csv(gff_file, sep="\t", header=None, comment='#')
    gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    # Filter the GFF file for features within the region of interest
    gff_filtered = gff_df[(gff_df['seqid'] == chromosome) & 
                          (gff_df['start'] >= start_pos) & 
                          (gff_df['end'] <= end_pos)].copy()

    # Adjust the start and end positions if start_pos > 1
    if start_pos > 1:
        gff_filtered['start'] = gff_filtered['start'] - (start_pos - 1)
        gff_filtered['end'] = gff_filtered['end'] - (start_pos - 1)

    # Save the filtered and adjusted GFF file as {Sample_Short}_{Chromosome}_extracted.gff
    extracted_gff = os.path.join(output_dir, f"{sample_short}_{chromosome}_extracted.gff")
    gff_filtered.to_csv(extracted_gff, sep="\t", header=False, index=False)
    print(f"Extracted GFF annotations written to {extracted_gff}")

    # Read the FASTA file and extract the relevant region
    output_fasta = os.path.join(output_dir, f"{sample_short}_{chromosome}_extracted.fasta")
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id == chromosome:
                # Extract the specific region and write it to a new FASTA file
                sub_record = record[start_pos-1:end_pos]  # Convert to 0-based indexing
                SeqIO.write(sub_record, output_fasta, "fasta")
    print(f"Extracted sequence written to {output_fasta}")
