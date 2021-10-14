#!/bin/bash 
#PBS -N snpEff
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=2:mem=16gb:scratch_local=8gb
#PBS -m abe
#PBS -j oe

module add jdk-8
#module add snpeff/snpeff-2017-11-24-intel-19.0.4-zxpnpan
wd=$PBS_O_WORKDIR
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
#cp -r /storage/plzen1/home/holcovam/programs/snpEff $SCRATCHDIR || exit 1
#cp $wd/$vcf $SCRATCHDIR || exit 1
#cd $SCRATCHDIR || exit 2
#echo data loaded at `date`

cp -r /storage/plzen1/home/holcovam/programs/snpEff $wd
cd $wd

java -Xmx16g -jar snpEff/snpEff.jar LyV2 ${vcf} > ${vcf/vcf.gz/ann.vcf}
rm $vcf
java -Xmx16g -jar snpEff/SnpSift.jar extractFields ${vcf/vcf.gz/ann.vcf} "CHROM" "POS" "GEN[*].GT" "ANN[*].GENEID" > ${vcf/vcf.gz/ann.tab.txt}

#rm -r snpEff
cp $SCRATCHDIR/* $wd/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
