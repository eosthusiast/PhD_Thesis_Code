#!/bin/bash -e
#SBATCH --job-name=herad_megahit_L01 # job name (shows up in the queue)
#SBATCH --time=10:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32G         # Memory in GB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config

cd /nesi/nobackup/ga03488/David_H/refreshed_backup/PPoutput/L01_first_run_20240212/trimmed/Ppo

# Load MEGAHIT module (adjust the module name as per your system)
module load MEGAHIT

# Define a function to run MEGAHIT on paired-end files
run_megahit2() {
    local read1=$1
    local read2=$2
    local output_prefix=$(basename ${read1} _1.fastq.gz)
    megahit -1 ${read1} -2 ${read2} -o ./megahit/${output_prefix}_megahit --num-cpu-threads 16 --memory 0.9
}

# Export the function to be used by parallel
export -f run_megahit2

# Find all paired-end files and run MEGAHIT on them
for read1 in *_1.fastq.gz; do
    read2=${read1/_1.fastq.gz/_2.fastq.gz}
    if [[ -f "$read2" ]]; then
        echo "Running MEGAHIT on $read1 and $read2"
        run_megahit2 "$read1" "$read2"
    else
        echo "Pair for $read1 not found, skipping."
    fi
done