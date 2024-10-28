#!/bin/bash

cd ./exonerate_2n_STE3_flank/

# Set the query file (CL4 and HD2 protein sequences) and the subject directory
CL4="../target_gene_regions/ICMP9630_FUN_003665_flank.fasta"
SUBJECT_DIR="/nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/05_FunFinder_Pangenome/input/*/"
ANNOTATED_GFF_DIR="/nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/06_HD1/protein_extraction/genomes_annotated/"

# Create a combined output text file
COMBINED_OUTPUT_FILE="combined_target_ids_and_columns.txt"
> $COMBINED_OUTPUT_FILE  # Clear the file before starting

# Loop through each genome in the subject directory
for GENOME in $SUBJECT_DIR/*_proteins.fasta; do
    # Extract the base name of the genome for output naming
    BASENAME=$(basename $GENOME _proteins.fasta)

    # Run exonerate using the p2g (protein-to-genome) model and capture all hits, filter out any line with '--' (like the exonerate completion line)
    exonerate -m affine:local -S FALSE -Q protein -T protein --query $CL4 --target $GENOME --showalignment no --showvulgar no --bestn 2 --ryo ">%ti\n%tas\n" | grep -v "^Command" | grep -v "^Hostname" | grep -v "^--" > "${BASENAME}_CL4_homologs_bestn.fasta"

    # Extract the full target sequences using the IDs from Exonerate results
    grep ">" "${BASENAME}_CL4_homologs_bestn.fasta" | sed 's/>//' | while read -r TARGET_ID; do
        # Extract the full protein sequence from the original target protein fasta
        awk -v id="$TARGET_ID" '/^>/{f=($0~id)}f' $GENOME >> "${BASENAME}_CL4_homologs_full.fasta"

        # Search for the target ID in the corresponding .gff file
        GFF_FILE="${ANNOTATED_GFF_DIR}${BASENAME}.gff"
        if [ -f "$GFF_FILE" ]; then
            # Extract the line containing $TARGET_ID and get the first column
            GFF_LINE=$(grep -m 1 "$TARGET_ID" "$GFF_FILE")
            if [ -n "$GFF_LINE" ]; then
                FIRST_COLUMN=$(echo "$GFF_LINE" | awk '{print $1}')
                
                # Write the first column and $TARGET_ID to the combined output file
                echo "$FIRST_COLUMN $TARGET_ID" >> $COMBINED_OUTPUT_FILE
            else
                echo "No match for $TARGET_ID in $GFF_FILE" >> error_log.txt
            fi
        else
            echo "GFF file $GFF_FILE not found" >> error_log.txt
        fi
    done

    echo "Finished processing $GENOME"
done

# Create a tar archive of the results
tar -cvf exo_flank.tar *
