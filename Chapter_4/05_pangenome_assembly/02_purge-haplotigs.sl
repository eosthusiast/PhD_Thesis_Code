#!/bin/bash -e
#SBATCH --job-name=02_haplotigs1 # job name (shows up in the queue)
#SBATCH --time=2:14:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

module purge
module load purge_haplotigs
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/01_purged_noragtag

# Directory paths
genome_dir="/nesi/nobackup/ga03488/David_H/Ppo_Popgen/03_ragtag_tuber-regium"

# Read through each fasta file in the genome_dir
for fasta_file in $genome_dir/*.fasta; do
  # Extract the basename (e.g., ICMP9630)
  basename=$(basename "$fasta_file" | cut -d'_' -f1)
  
  # Extract the l, m, h values from cov_cutoffs_all28.txt
  read l m h <<< $(awk -v id="$basename" '$1 == id {print $2, $3, $4}' ../cov_cutoffs_all28.txt)
  
  # Ensure the values were found
  if [ -z "$l" ] || [ -z "$m" ] || [ -z "$h" ]; then
    echo "No coverage cutoffs found for $basename. Skipping."
    continue
  fi
  
  # Find the corresponding BAM file
  bam_file="${basename}_selfmap.bam.gencov"
  echo "Processing basename: $basename"
  echo "Expected BAM file: $bam_file"
  echo "l: $l, m: $m, h: $h"
  
  if [ -f "$bam_file" ]; then
    # Define the output file name for coverage stats
    coverage_stats="${basename}_coverage_stats.csv"
    
    # Run purge_haplotigs cov with the correct parameters and output file
    echo "Running purge_haplotigs cov for $basename with l=$l, m=$m, h=$h"
    purge_haplotigs cov -i "$bam_file" -l "$l" -m "$m" -h "$h" -o "$coverage_stats"
    
    # Run purge_haplotigs purge with the relevant genome and coverage stats file
    echo "Running purge_haplotigs purge for $basename"
    purge_haplotigs purge -g "$fasta_file" -c "$coverage_stats" -t 16 -o "${basename}"
  else
    echo "BAM file $bam_file not found for $basename. Skipping."
  fi
done


# Generate dotplots from final alignments - option: -a / -align_cov     Percent cutoff for identifying a contig as a haplotig. DEFAULT = 70
#purge_haplotigs purge  -g ../trialDH641/annotate_results/Pleurotus_purpureo-olivaceus_DH641.scaffolds.fa  -c coverage_stats.csv  -a 70  -d  -b PH.bam

