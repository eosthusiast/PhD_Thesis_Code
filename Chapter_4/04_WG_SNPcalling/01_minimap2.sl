#!/bin/bash -e
#SBATCH --job-name=02_mapping # job name (shows up in the queue)
#SBATCH --time=10:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/14_whole_genome_snpcalling/

# Load necessary modules
module purge
module load minimap2
module load SAMtools

ref=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/haploidDH641_noragtag/annotate_results/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa #Reference genome for alignment
datadir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/01_purged_noragtag/curated_genomes/done/ #Directory with fastq data
samdir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/14_whole_genome_snpcalling/01_mapping/sam/ #Sam file output
bamdir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/14_whole_genome_snpcalling/01_mapping/bam/ #Bam file output
vardir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/14_whole_genome_snpcalling/02_variants/

# Loop through all .fasta files in the datadir
#for fasta in "$datadir"/*.fasta; do
#    basename=$(basename "$fasta" .fasta)
#
#   # Align fasta to reference genome using minimap2 and output SAM file
#    minimap2 -ax asm5 $ref "$fasta" > "${samdir}/${basename}.sam"
#
#    # Convert SAM to BAM, sort, and index
#    samtools view -bS "${samdir}/${basename}.sam" | samtools sort -o "${bamdir}/${basename}.sorted.bam"
#    samtools index "${bamdir}/${basename}.sorted.bam"
#
#    # Clean up SAM file
#    rm "${samdir}/${basename}.sam"
#done

# Create a list of BAM files
#ls $bamdir/*.sorted.bam > $bamdir/bam_list.txt

module load BCFtools

# Set output VCF file
OUTPUT_VCF=$vardir/raw_variants.vcf
RENAMED_VCF=$vardir/raw_variants_renamed.vcf

# Generate pileup and call variants
bcftools mpileup -a AD,DP,SP -f $ref -b $bamdir/bam_list.txt -Ob -o $vardir/all_samples_annotated.bcf --threads 16
bcftools call -f GQ,GP -m -v -Ov -o $OUTPUT_VCF $vardir/all_samples_annotated.bcf

### Rename sample names with sample labels, rather than full bam filepaths and extensions

# Check that both files exist
if [ ! -f $bamdir/bam_list.txt ] || [ ! -f $vardir/new_names.txt ]; then
    echo "One or both of the files sample_names.txt or new_names.txt do not exist."
    exit 1
fi

# Copy the original VCF to the new output VCF
cp $OUTPUT_VCF $RENAMED_VCF

# Read both files line by line and replace the old sample names with the new ones in the output VCF
paste $bamdir/bam_list.txt $vardir/new_names.txt | while IFS=$'\t' read -r old_name new_name; do
    sed -i "s|$old_name|$new_name|g" $RENAMED_VCF
done

echo "Sample names have been replaced and saved in $RENAMED_VCF."