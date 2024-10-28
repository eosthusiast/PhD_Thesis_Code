from Bio import SeqIO
import pandas as pd
import os
import glob

# Define the input and output directories
genome_dir = "./genomes_annotated"
output_dir = "./STE3_flank_ICMP9630_n2"
os.makedirs(output_dir, exist_ok=True)

# Load the matB_gff_table.txt file
matB_gff_table = pd.read_csv("STE3flank_target_ids_and_columns-Copy1.txt", sep="\t")

# Print the first few rows and columns to verify loading is correct
print(matB_gff_table.head())
print(matB_gff_table.columns)

# Loop through each row in the matB_gff_table
for idx, row in matB_gff_table.iterrows():
    # Extract the sample from the 'Chromosome' column
    sample = row['Chromosome'].split('_')[0]
    protein_name = row['ProteinName']
    chromosome = row['Chromosome']

    # Use wildcards to search for the FASTA and GFF files
    fasta_file = glob.glob(os.path.join(genome_dir, f"*{sample}*.fa"))
    gff_file = glob.glob(os.path.join(genome_dir, f"*{sample}*.gff3"))

    # Check if the files were found
    if len(fasta_file) == 0:
        print(f"FASTA file for sample {sample} not found, skipping...")
        continue
    if len(gff_file) == 0:
        print(f"GFF file for sample {sample} not found, skipping...")
        continue

    # Use the first match found for both FASTA and GFF
    fasta_file = fasta_file[0]
    gff_file = gff_file[0]

    # Create output files
    output_fasta = os.path.join(output_dir, f"{sample}_{chromosome}_flank.fasta")
    output_gff = os.path.join(output_dir, f"{sample}_{chromosome}_flank.gff")

    # Extract the entire chromosome from the FASTA file
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id == chromosome:
                # Set start as 1 and end as the length of the chromosome
                start = 1
                end = len(record)
                
                # Extract the chromosome sequence (entire sequence)
                sub_record = record[start-1:end]  # 0-based indexing in Python
                
                # Write the modified sequence to the output FASTA
                SeqIO.write(sub_record, output_fasta, "fasta")
                print(f"Extracted chromosome {chromosome} for {protein_name} into {output_fasta}")

    # Filter the GFF file to extract annotations for the entire chromosome
    gff_df = pd.read_csv(gff_file, sep="\t", header=None, comment='#')
    gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    # Filter GFF entries for the entire chromosome
    gff_filtered = gff_df[gff_df['seqid'] == chromosome]

    # Write the filtered GFF to a new file
    gff_filtered.to_csv(output_gff, sep="\t", header=False, index=False)
    print(f"Extracted GFF for {chromosome} in {protein_name} into {output_gff}")
