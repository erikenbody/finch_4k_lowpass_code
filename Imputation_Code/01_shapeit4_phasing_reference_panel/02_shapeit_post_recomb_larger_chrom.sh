#!/bin/bash
#SBATCH -A snic2021-22-187
#SBATCH -p node
#SBATCH -C mem256GB #1TB was needed for chr2
#SBATCH -M rackham #rackham only
#SBATCH -t 5-00:00:00 #should have run for longer than 5 days (1 and 2 finished though)
#SBATCH -J shapeit4
#SBATCH -e shapeit4_%J_%A_%a.err
#SBATCH -o shapeit4_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

##for chromosome 2 this needed the 1TB node previously. So lets just run them all on that big node
##seems to max out at about 10 threads
##ran tests just with a small chromosome and running this as an array
#original notes were: #6 for 1-31, 256 needed for first 10 chr or so. 1TB for array 2

module load bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.12

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/23_shapeit4_phasing_new_samps
cd $WORK_D

##using 4.2.1 - should be faster. Install required editing the makefile as before
SI=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/shapeit4-4.2.1/bin
INTERVAL_LIST=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt
SAMPLE_DIR=$WORK_D

#for SLURM_ARRAY_TASK_ID in {1..3}
for SLURM_ARRAY_TASK_ID in 3
do
  INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
  #gather all whatshap vcfs (two rounds, one before new 2020 samples and one after)
  ls -1 $SAMPLE_DIR/*/*${INTERVAL}.whatshap.vcf.gz > ${INTERVAL}_415whatshap.list
  bcftools merge -l ${INTERVAL}_415whatshap.list -O b -o ${INTERVAL}.415samps.bcf --threads 20
  bcftools index ${INTERVAL}.415samps.bcf --threads 20

  $SI/shapeit4.2 --input ${INTERVAL}.415samps.bcf \
            --output ${INTERVAL}.415samps.postrecomb.bcf \
            --sequencing \
            --map gmap_files/${INTERVAL}_parv.gmap \
            --thread 20 \
            --region $INTERVAL \
            --log ${INTERVAL}.phasing.postrecomb.log

  bcftools index ${INTERVAL}.415samps.postrecomb.bcf
done

# for SLURM_ARRAY_TASK_ID in {1..31}
# do
#   INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
#   ls -1 ${INTERVAL}.415samps.postrecomb.bcf
# done
