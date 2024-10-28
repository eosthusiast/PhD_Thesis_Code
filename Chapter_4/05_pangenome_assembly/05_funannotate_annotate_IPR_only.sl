#!/bin/bash -e
#SBATCH --job-name=13_iprsc_funannotate # job name (shows up in the queue)
#SBATCH --time=11:13:13      # Walltime (HH:MM:SS)
#SBATCH --mem=150GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49254984 # wait till previous task is complete
#SBATCH --mail-user=david.hera@pg.canterbury.ac.nz
#SBATCH --mail-type=ALL

module purge
module load Apptainer
module load GeneMark-ES
module load InterProScan/5.51-85.0-gimkl-2020a-Perl-5.30.1-Python-3.8.2
module load PCRE2
#module load eggnog-mapper
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/12_pangenome_prep_all_Ppo/03_annotated_genomes

# Loop over each .fasta file in the curated_genomes directory
for fasta_file in ../01_purged_noragtag/curated_genomes/13/*.fasta; do
    # Extract the basename (without extension)
    basename=$(basename "$fasta_file" .fasta)
    
    # Define the input and output files
    input_file="$fasta_file"
    output_file="../01_purged_noragtag/curated_genomes/${basename}_renamed.fasta"
    
    if [ ! -d "${basename}_haploid/annotate_misc/" ]; then
    mkdir -p "${basename}_haploid/annotate_misc/"
    fi

    /opt/nesi/CS400_centos7_bdw/InterProScan/5.51-85.0-gimkl-2020a-Perl-5.30.1-Python-3.8.2/interproscan.sh -i ./"${basename}_haploid"/predict_results/*.proteins.fa --cpu 15 --output-dir "${basename}_haploid"/annotate_misc/ -goterms -pa -f XML 2>&1 | tee -a "iprscan3-${basename}_haploid.log"

    # rename interproscan results
    mv ./"${basename}_haploid"/annotate_misc/*.xml ./"${basename}_haploid"/annotate_misc/iprscan.xml
    
    echo "IPRScan finished on ${basename}, moving on to next sample now."
    
    #apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate remote -i "${basename}_haploid" -m antismash -e david.hera@pg.canterbury.ac.nz

    # Run EggNog and Phobius manually on webserver, add to results folder, then run next code
    #echo "Run EggNog manually on webserver, add to results folder, then run next code"

    #apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate annotate -i haploidDH641_noragtag --cpus 16 --phobius ./haploidDH641_noragtag/annotate_misc/phobius_result.txt

done


