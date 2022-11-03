#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 10-00:00:00
#SBATCH -J map_lowpass
#SBATCH -e map_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o map_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools samtools bwa/0.7.17

#run on chunks of 100 samples
#20 arrays

REF_DIR=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA
REF=Camarhynchus_parvulus_V1.0.fasta
SAMPLE_DIR=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/lowpass/RAW_FASTQ/novaseq_run2/SH-2268/190906_A00181_0112_BHLGTVDSXX
WORK_D=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/lowpass/lowpass_pipeline/01.2_mapping

#do this once:
# cd $SAMPLE_DIR
# ls -d -1 $PWD/*/ > $WORK_D/INPUT_DIRS.txt
# cd $WORK_D
# split -l 30 -d INPUT_DIRS.txt finch_split.

cd $WORK_D

SAMPLES=`ls -1 *split* | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

cp $REF_DIR/Camarhynchus_parvulus_V1.0.fasta $SNIC_TMP
cp $REF_DIR/Camarhynchus_parvulus_V1.0.dict $SNIC_TMP
cp $REF_DIR/Camarhynchus_parvulus_V1.0.fasta.* $SNIC_TMP

while read FILENAME; do
  SAMPLE1=$(basename $FILENAME)
  SAMPLE=${SAMPLE1/Sample_SH-2268-/}
  R1=${FILENAME}*R1*fastq.gz
  R2=${FILENAME}*R2*fastq.gz
  bR1=$(basename $R1)
  bR2=$(basename $R2)
  cp $R1 $SNIC_TMP
  cp $R2 $SNIC_TMP
  if [ -f ${SAMPLE}.sort.bam ]; then echo "bam exists"; else cp $R1 $SNIC_TMP ; cp $R2 $SNIC_TMP ; fi
  if [ -f ${SAMPLE}.sort.bam ]; then echo "bam exists"; else bwa mem -t 8 -R $(echo "@RG\\tID:$SAMPLE\\tSM:$SAMPLE") $SNIC_TMP/$REF $SNIC_TMP/$bR1 $SNIC_TMP/$bR2 | samtools sort -@ 8 -O bam -o $SNIC_TMP/${SAMPLE}.sort.bam ; samtools index -@ 8 $SNIC_TMP/${SAMPLE}.sort.bam ; mv $SNIC_TMP/*bam* .; fi
done < $SAMPLES
