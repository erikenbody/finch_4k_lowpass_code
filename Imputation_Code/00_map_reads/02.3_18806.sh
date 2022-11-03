#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 00-05:00:00
#SBATCH -J map_lowpass
#SBATCH -e map_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o map_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#this sample didnt get mapped

module load bioinfo-tools samtools bwa/0.7.17

WORK_D=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/lowpass/lowpass_pipeline/01.2_mapping
cd $WORK_D

bR1=/crex/proj/uppstore2017190/private/02_RAW_RESEQUENCING_DATA/finches_lowpass_raw_data_2019/RAW_FASTQ/novaseq_run2/SH-2268/190906_A00181_0112_BHLGTVDSXX/Sample_SH-2268-01Dap18806/SH-2268-01Dap18806_S239_L001_R1_001.fastq.gz
bR2=/crex/proj/uppstore2017190/private/02_RAW_RESEQUENCING_DATA/finches_lowpass_raw_data_2019/RAW_FASTQ/novaseq_run2/SH-2268/190906_A00181_0112_BHLGTVDSXX/Sample_SH-2268-01Dap18806/SH-2268-01Dap18806_S239_L001_R2_001.fastq.gz

REF=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta

bwa mem -t 8 -R $(echo "@RG\\tID:01Dap18806\\tSM:01Dap18806") $REF $bR1 $bR2 | samtools sort -@ 8 -O bam -o 01Dap18806.sort.bam ; samtools index -@ 8 01Dap18806.sort.bam
