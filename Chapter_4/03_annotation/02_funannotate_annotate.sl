#!/bin/bash -e
#SBATCH --job-name=herad_funannotate # job name (shows up in the queue)
#SBATCH --time=10:23:13      # Walltime (HH:MM:SS)
#SBATCH --mem=32GB          # Memory in MB
# #SBATCH --qos=debug          # debug QOS for high priority job tests
#SBATCH --account=ga03488 # define project code
#SBATCH --cpus-per-task=16 # number of logical CPUs per task, matching the n THREADZ in pp_align.config
# #SBATCH --dependency=afterok:49088000 # wait till previous purge haplotype task is complete

module purge
module load Apptainer
module load GeneMark-ES
#module load InterProScan/5.51-85.0-gimkl-2020a-Perl-5.30.1-Python-3.8.2
module load PCRE2
#module load eggnog-mapper
cd /nesi/nobackup/ga03488/David_H/Ppo_Popgen/11_funannotate_predict

#/opt/nesi/CS400_centos7_bdw/InterProScan/5.51-85.0-gimkl-2020a-Perl-5.30.1-Python-3.8.2/interproscan.sh -i ./haploidDH641/predict_results/*.proteins.fa --cpu 15 --output-dir ./haploidDH641/annotate_misc/ -goterms -pa -f XML 2>&1 | tee -a iprscan3-bigmem_hap.log

# rename interproscan results
#mv ./haploidDH641/annotate_misc/*.xml ./haploidDH641/annotate_misc/iprscan.xml

#apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate remote -i haploidDH641 -m antismash -e david.hera@pg.canterbury.ac.nz

# Run EggNog manually on webserver, add to results folder, then run next code

apptainer run /nesi/nobackup/ga03488/David_H/funannotate/funannotate.sif funannotate annotate -i haploidDH641 --cpus 16


