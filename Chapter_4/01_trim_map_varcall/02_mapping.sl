#!/bin/bash -e
#SBATCH --job-name=02_mapping # job name (shows up in the queue)
#SBATCH --time=20:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=150GB         # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

# Load necessary modules
module purge

module load BWA
module load SAMtools
module load BCFtools

ref=/nesi/nobackup/ga03488/David_H/Ppo_Popgen/13_DH641_haploid_non_ragtag_assembly/haploidDH641_noragtag/annotate_results/Pleurotus_purpureo-olivaceus_DH641_haploid_noragtag.scaffolds.fa #Reference genome for alignment
datadir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/01_trimmed/ #Directory with fastq data
samdir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/02_mapped/sam/ #Sam file output
bamdir=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/02_mapped/bam/ #Bam file output
fq1=_1_val_1.fq.gz #Read 1 suffix
fq2=_2_val_2.fq.gz #Read 2 suffix

#First index the reference genome
#time bwa index $ref

#Now, map to reference genome
#for samp in ${datadir}*_1_val_1.fq.gz #Remember to be explicit with file location
#do
#    base=$(basename ${samp} _1_val_1.fq.gz)
#    
#    echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
#    time bwa mem -M -t 16 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam
#
#    echo "Converting sam file to bam file for $base"
#    time samtools view -T $ref -b ${samdir}${base}.sam > ${bamdir}${base}.bam
#done

# 

data=/nesi/nobackup/ga03488/David_H/Ppo_Popgen_retrimmed/02_mapped/bam/

# Marking duplicates (and removing them with -r)
#for bam in ${data}/*.bam
#do
#    base=$(basename ${bam} .bam)
#    echo "Now preparing to mark duplicates for ${base}..."
#    samtools sort -@ 8 -n -o ${data}nodup_bam/${base}.nsorted.bam ${bam}
#    samtools fixmate -@ 8 -r -m -c ${data}nodup_bam/${base}.nsorted.bam \
#        ${data}nodup_bam/${base}.fixmate.bam
#    samtools sort -@ 8 -o ${data}nodup_bam/${base}.fixmate.sorted.bam \
#        ${data}nodup_bam/${base}.fixmate.bam
#    samtools markdup -r -@ 8 ${data}nodup_bam/${base}.fixmate.sorted.bam \
#        ${data}nodup_bam/${base}_nodup.bam
#    samtools stats ${bam} > ${data}nodup_bam_stats/${base}.stats
#    samtools index ${data}nodup_bam/${base}_nodup.bam
#done

module load mosdepth

# quality metrics with mosdepth - removed qualimap from Jana's script because it's not installed

#for bam in ${data}nodup_bam/*_nodup.bam
#    do
#    base=$(basename ${bam} _nodup.bam)
#    echo "Running Qualimap for ${base}..."
#    echo "Running calculating stats for ${base}..."
#    mosdepth --threads 16 --fast-mode --by 50 ${data}nodup_bam_stats/${base} ${bam}
#done

# mosdepth plot
python plot-dist.py ${data}nodup_bam_stats/*.global.dist.txt