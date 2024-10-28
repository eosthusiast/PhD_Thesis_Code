#!/bin/bash

cd ./exonerate_2n_fulllengthprotein/

# Set the query file (HD1 and HD2 protein sequences) and the subject directory
HD1="../target_gene_regions/DH641_k141_8609_HD1_protein.fasta"
HD2="../target_gene_regions/DH641_k141_8609_HD2_protein.fasta"
HD2short="../target_gene_regions/DH639_k141_HD2_short_protein.fasta"
STE3_1="../target_gene_regions/DH641_k141_STE3_1.fasta"
SUBJECT_DIR="/nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/05_FunFinder_Pangenome/input/*/"

# Loop through each genome in the subject directory
for GENOME in $SUBJECT_DIR/*_proteins.fasta; do
    # Extract the base name of the genome for output naming
    BASENAME=$(basename $GENOME _proteins.fasta)

    # Run exonerate using a affine:local p2p (protein-to-protein) model and capture all hits, filter out any line with '--' (like the exonerate completion line)
    exonerate -m affine:local -S FALSE -Q protein -T protein --query $HD1 --target $GENOME --showalignment no --showvulgar no --bestn 2 --ryo ">%ti\n%tas\n" | grep -v "^Command" | grep -v "^Hostname" | grep -v "^--" > "${BASENAME}_HD1_homologs_bestn.fasta" 

    # Extract the full target sequences using the IDs from Exonerate results
    grep ">" "${BASENAME}_HD1_homologs_bestn.fasta" | sed 's/>//' | while read -r TARGET_ID; do
        # Extract the full protein sequence from the original target protein fasta
        awk -v id="$TARGET_ID" '/^>/{f=($0~id)}f' $GENOME >> "${BASENAME}_HD1_homologs_full.fasta"
    done

    exonerate -m affine:local -S FALSE -Q protein -T protein --query $HD2 --target $GENOME --showalignment no --showvulgar no --bestn 2 --ryo ">%ti\n%tas\n" | grep -v "^Command" | grep -v "^Hostname" | grep -v "^--" > "${BASENAME}_HD2_homologs_bestn.fasta" 
    
    # Extract full sequences for HD2
    grep ">" "${BASENAME}_HD2_homologs_bestn.fasta" | sed 's/>//' | while read -r TARGET_ID; do
        awk -v id="$TARGET_ID" '/^>/{f=($0~id)}f' $GENOME >> "${BASENAME}_HD2_homologs_full.fasta"
    done

    exonerate -m affine:local -S FALSE -Q protein -T protein --query $HD2short --target $GENOME --showalignment no --showvulgar no --bestn 2 --ryo ">%ti\n%tas\n" | grep -v "^Command" | grep -v "^Hostname" | grep -v "^--" > "${BASENAME}_HD2short_homologs_bestn.fasta" 

    # Extract full sequences for HD2short
    grep ">" "${BASENAME}_HD2short_homologs_bestn.fasta" | sed 's/>//' | while read -r TARGET_ID; do
        awk -v id="$TARGET_ID" '/^>/{f=($0~id)}f' $GENOME >> "${BASENAME}_HD2short_homologs_full.fasta"
    done
    
    exonerate -m affine:local -S FALSE -Q protein -T protein --query $STE3_1 --target $GENOME --showalignment no --showvulgar no --bestn 2 --ryo ">%ti\n%tas\n" | grep -v "^Command" | grep -v "^Hostname" | grep -v "^--" > "${BASENAME}_STE3_1_homologs_bestn.fasta" 

    # Extract full sequences for HD2short
    grep ">" "${BASENAME}_STE3_1_homologs_bestn.fasta" | sed 's/>//' | while read -r TARGET_ID; do
        awk -v id="$TARGET_ID" '/^>/{f=($0~id)}f' $GENOME >> "${BASENAME}_STE3_1_homologs_full.fasta"
    done

    echo "Finished processing $GENOME"
    
    tar -cvf exoall.tar *
done
