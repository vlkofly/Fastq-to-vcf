#!/bin/bash
#PBS -N merge_vcf
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=4:mem=50gb:scratch_local=100gb
#PBS -j oe

# this is a continuation of filtering F, actually merging was done wrong there so when this script will be finished, 
# I can only merge the good and part of F filtering with this script.

if [ -z "$vcf" ]; then
        echo "Error! Specify a key string common for all vcfs you want to merge "
        echo "qsub -v 'vcf=blablabla.vcf.gz'"
        exit 1
        fi

if [ -z "$out" ]; then
        echo "Specify name of the outfile"
	echo " if not specified default will be used"
	out=${vcf}.merged.vcf.gz
        fi

if [ -z "$ref" ]; then
        echo "Specify reference full path, default reference used instead"
        echo "qsub -v '/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta'"
        ref=/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta
        fi
wd=$PBS_O_WORKDIR
nt=4 # number of threads change based on number of scaffolds - but it is better to run one job per a scaffold

refdir=`dirname $ref`
ref=`echo $ref | sed 's/^.*\/\(.*\/.*\)$/\1/g'`

trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 


cp -r $refdir $SCRATCHDIR/
cp $wd/${vcf} $SCRATCHDIR/ # modified 12.2.2021 to merge pops
#cp $wd/*${vcf}*vcf.gz $SCRATCHDIR/
#cp $wd/*${vcf}*vcf.gz.tbi $SCRATCHDIR/ # copy index as well

cd $SCRATCHDIR
ls

module add gatk-3.7
module add bcftools-1.6
module add htslib-1.9
tmj="-Djava.io.tmpdir=$SCRATCHDIR"

# in the case of vcf originating from Scantools you must zip it
ls *vcf | parallel -j $nt "bgzip -c {} > {}.gz; tabix -p vcf {}.gz"


mergelist=""

for f in *.vcf.gz
do
mergelist+="-V $f "
done

echo $mergelist

java $tmj -cp $GATK/GenomeAnalysisTK.jar \
	org.broadinstitute.gatk.tools.CatVariants -R $ref \
	$mergelist -assumeSorted -out $out

####stats
echo "vcf_file numvariants"
for v in `ls *vcf*gz`
do
n=`zcat $v | grep -v "#" | wc -l`
echo $v $n
done >> ${out/vcf.gz/numbers.ssv}

cat ${out/vcf.gz/numbers.ssv}

  
cp $SCRATCHDIR/$out ${wd}/
cp ${out/vcf.gz/numbers.ssv} ${wd}/ || export CLEAN_SCRATCH=false

echo "filtering of ${vcf} done"

