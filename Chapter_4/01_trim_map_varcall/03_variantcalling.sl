#!/bin/bash -e
#SBATCH --job-name=03_varcall # job name (shows up in the queue)
#SBATCH --time=3:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

REF_GENOME=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/haploidDH641_noragtag/annotate_results/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa
TRIMMED_FASTQ_DIR=/nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L01_first_run_20240212/trimmed/Ppo/trimmed_fastq
OUTPUT_DIR=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/02_mapped/bam/nodup_bam
VAR_DIR=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/03_variant_calling

cd $OUTPUT_DIR

# Load necessary modules
module purge
module load BCFtools

# Create a list of BAM files
ls $OUTPUT_DIR/*_nodup.bam > $OUTPUT_DIR/bam_list.txt

# Create an output VCF file
OUTPUT_VCF=$VAR_DIR/raw_variants.vcf
RENAMED_VCF=$VAR_DIR/raw_variants_renamed.vcf

# Generate pileup and call variants
#bcftools mpileup -a AD,DP,SP -f $REF_GENOME -b bam_list.txt -Ob -o all_samples.bcf --threads 16
bcftools call -f GQ,GP -m -v -Ov -o $OUTPUT_VCF all_samples.bcf --threads 16

### Rename sample names with sample labels, rather than full bam filepaths and extensions

# Check that both files exist
if [ ! -f bam_list.txt ] || [ ! -f $VAR_DIR/new_names.txt ]; then
    echo "One or both of the files sample_names.txt or new_names.txt do not exist."
    exit 1
fi

# Copy the original VCF to the new output VCF
cp $OUTPUT_VCF $RENAMED_VCF

# Read both files line by line and replace the old sample names with the new ones in the output VCF
paste bam_list.txt $VAR_DIR/new_names.txt | while IFS=$'\t' read -r old_name new_name; do
    sed -i "s|$old_name|$new_name|g" $RENAMED_VCF
done

echo "Sample names have been replaced and saved in $RENAMED_VCF."

