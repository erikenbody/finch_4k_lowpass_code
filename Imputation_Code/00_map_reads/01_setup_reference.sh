#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -p core -n 4
#SBATCH -t 0-23:00:00
#SBATCH -J dict
#SBATCH -e dict_%A_%a.err            # File to which STDERR will be written
#SBATCH -o dict_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com


module load bwa/0.7.17 picard/2.10.3

cd /crex/proj/uppstore2017190/private/01_REFERENCE_DATA

REF=Camarhynchus_parvulus_V1.0.fasta

bwa index $REF

java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary \
      R=$REF \
      O=Camarhynchus_parvulus_V1.0.dict
