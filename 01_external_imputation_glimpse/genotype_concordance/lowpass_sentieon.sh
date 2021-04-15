#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
data_dir=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/test_sample

# Update with the location of the reference data files
fasta=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta

#testing install in new place
# Set SENTIEON_LICENSE if it is not set in the environment
#export SENTIEON_LICENSE=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools/sentieon/Uppsala_University_eval.lic
export SENTIEON_LICENSE=/crex/proj/snic2020-2-19/private/wagtail/tools/sentieon/Uppsala_University_eval.lic
# Update with the location of the Sentieon software package
#SENTIEON_INSTALL_DIR=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools/sentieon/sentieon-genomics-201911
SENTIEON_INSTALL_DIR=/home/eenbody/snic2020-2-19/private/wagtail/tools/sentieon/sentieon-genomics-201911
# Update with the location of temporary fast storage and uncomment

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

sample=$1
bam_input=$2

# Other settings
nt=5 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$data_dir/$sample #ede edited this to be sample
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 6. HC Variant caller
# Note: Sentieon default setting matches versions before GATK 3.7.
# Starting GATK v3.7, the default settings have been updated multiple times.
# Below shows commands to match GATK v3.7 - 4.1
# Please change according to your desired behavior.
# ******************************************



# Matching GATK 4.1
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $bam_input --algo Haplotyper --genotype_model multinomial --emit_mode all --emit_conf 30 --call_conf 30 ${sample}.g.vcf.gz

#remove confidence
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $bam_input --algo Haplotyper --genotype_model multinomial --emit_mode confident --emit_conf 0 --call_conf 0 ${sample}.vcf.gz
