#!/bin/bash
#PBS -N index_ref
#PBS -l select=1:ncpus=1:mem=50gb:scratch_local=50gb
#PBS -l walltime=8:00:00
#PBS -j oe

# this script takes reference genome fasta file and prepares
# it for aligning with bwa, indices and dictionaries are built by several
# programs. 
# input reference fasta should be specified by variable ref
trap 'clean_scratch' TERM EXIT

wd=$PBS_O_WORKDIR
#ref="9172_ref_Sturnus_vulgaris-1.0_chrUn.fa.gz" # user input
outdir=bwa
cp ${wd}/${ref} $SCRATCHDIR
cd $SCRATCHDIR
ls

module add picard-2.8.1
module add samtools-1.4
module add bwa-0.7.15
module add htslib-1.6
##bgzip deprecated

#if [[ $ref != *.gz ]]
#then
#bgzip -c $ref > ${ref}.gz
#ref=${ref}.gz
#fi


java -jar $PICARD281 CreateSequenceDictionary R=${ref} O=${ref/fasta/dict}

samtools faidx ${ref}

bwa index ${ref}

mkdir $outdir
mv * $outdir/

cp -r $outdir ${wd}/ || export CLEAN_SCRATCH=false

