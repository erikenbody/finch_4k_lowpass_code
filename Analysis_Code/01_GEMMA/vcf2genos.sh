#!/bin/bash
module load bcftools/1.12

VCF=${1}
SAMPLE=${2}
bcftools query -l $VCF | tr "\n" "\t" > ${SAMPLE}_header.txt
echo "" >> ${SAMPLE}_header.txt #add new line at end of header
bcftools query -f '[%GT\t]\n' $VCF > ${SAMPLE}_raw.genos
cat ${SAMPLE}_header.txt ${SAMPLE}_raw.genos > ${SAMPLE}.genos
bcftools query -f '%CHROM\t %POS\t %REF\t %ALT\n' $VCF > ${SAMPLE}_raw.pos
echo -e "CHROM\tPOS\tREF\tALT" | cat - ${SAMPLE}_raw.pos > ${SAMPLE}.pos
paste -d "\t" ${SAMPLE}.pos ${SAMPLE}.genos > ${SAMPLE}.pos.genos

#on mac:
#sed -i '' -e "s/0|1/1/g" -e "s/1|0/1/g" -e "s/0|0/0/g" -e "s/1|1/2/g" ${SAMPLE}.pos.genos
sed ${SAMPLE}.pos.genos -e "s:0|1:1:g" -e "s:1|0:1:g" -e "s:0|0:0:g" -e "s:1|1:2:g" > ${SAMPLE}.pos.bin.genos

rm ${SAMPLE}.pos
rm ${SAMPLE}_raw.pos
rm ${SAMPLE}.genos
rm ${SAMPLE}_raw.genos
rm ${SAMPLE}_header.txt
