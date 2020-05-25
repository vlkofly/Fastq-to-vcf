#!/bin/bash
#PBS -N V3_join_genotype
#PBS -l walltime=256:00:00
#PBS -l select=1:ncpus=1:mem=850gb:scratch_local=1200gb
#PBS -j oe

# modify number of cpus and other parameters according to your data
# use the bash script instead of the python as this is much shorter
# only you have to supply proper variables in the scripts, folders with data
# reference folder and output folder

# parallel version that splits genotypeGVCF per scaffold and then merges the vcfs together
# scatter-gather fix of parallelization
# 7.3.2019 script modified for more general use, where individual parameters will be 

### if clause for checking the input parameters TODO

if [ -z "$samples" ]; then
        echo "Error! Specify sample list (relative path from working directory) one column with new sample names!"
        echo "qsub -v 'samples=blablabla.txt'"
        exit 1
        fi

if [ -z "$outvcf" ]; then
        echo "specify output name, default name used instead"
        echo "qsub -v 'outvcf=blablabla.vcf.gz'"
        outvcf=${samples}.vcf.gz
        fi


if [ -z "$outdir" ]; then
        echo "specify output name, default name used instead"
        echo "qsub -v 'outvcf=blablabla.vcf.gz'"
        outdir="vcf_final"
        fi

if [ -z "$sourcedirs" ]; then
        echo "specify where g.vcf.gz files are, fullpath, different directories separated by space"
	echo "the sourcedirs has to have the permissions change with chgrp and chmod"
        echo "qsub -v 'sourcedirs=\"/a/b/c /d/g\" ' "
        exit 1
        fi

if [ "$NV" == "yes" ]; then # if you suply some text in variable NV, gatk will do calling of all sites not only variants
	echo "calling all sites"
	NV="--includeNonVariantSites"
else
	NV=""

	fi

if [ -z "$ref" ]; then
	echo "specify reference full path, default reference used instead"
	echo "qsub -v '/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta'"
	ref=/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta
	fi


if [ -z "$sc" ]; then
        echo "specify full path to scaffold list to run across in parallel"
	echo "if not specified default for lyrata is used"
        echo "one scaffold/chromosome per line"
        sc=/storage/brno3-cerit/home/filip_kolar/lyrata_parallel_scaffold_list.txt
        fi




#outvcf=halleri_snp_raw_repaired.vcf.gz
# the sourcedirs has to have the permissions change with chgrp and chmod
#sourcedirs="/storage/brno3-cerit/home/filip_kolar/HCrun2 /storage/brno3-cerit/home/vlkofly/HC300 /storage/brno3-cerit/home/vlkofly/HCmajda /storage/brno7-cerit/home/vlkofly/HCrun1 /storage/plzen1/home/vlkofly/HCrun0new /storage/brno3-cerit/home/veronikakonecna/HCserp" 
wdir=$PBS_O_WORKDIR
refdir=`dirname $ref`
ref=`echo $ref | sed 's/^.*\/\(.*\/.*\)$/\1/g'`
scaf=`cat $sc`
#samples="gvcf_halleri_samples.txt"
nt=1

module add gatk-3.7
module add parallel-20160622

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 

if [ ! -d "${wdir}/${outdir}" ] ; then echo "output dir created";  mkdir ${wdir}/${outdir}; fi # create output directory if not there yet

cp -r $refdir $SCRATCHDIR/

cp $sc $SCRATCHDIR # supply scaffold information in this file
cp /storage/brno3-cerit/home/filip_kolar/lyrata_scaffolds_rest.intervals $SCRATCHDIR # if the list is too long merge shorter scaffolds and upload separately


supplied_num=`wc -l ${wdir}/${samples}`

while read s # copy samples from various directories to scrach
do
find $sourcedirs -name ${s}* -type f -maxdepth 1 -exec cp {} $SCRATCHDIR \;
done < ${wdir}/${samples} 

echo "copying done `date`"


cd $SCRATCHDIR
ls
# here I will have to do some sample selection according the project

vf=`ls *g.vcf.gz`
samplenum=`echo $vf | tr " " \\n | wc -l`
#echo $vf

echo "genotyping $samplenum samples from $supplied_num samples supplied"

for f in $vf
do
gvlist+=" -V $f"
done 
#echo $gvlist
tmj="-Djava.io.tmpdir=$SCRATCHDIR"
#genotype GVCF
# consider adding -Xmx50g or some other value, if it fails again
java -XX:ParallelGCThreads=1 -Xmx45g -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R $ref $gvlist --max_num_PL_values 1500 $NV -nt 1 -L $scaf -o ${scaf}$outvcf 

rm  *g.vcf*

#echo "genotypeGVCF done in parallel for separate scaffolds:"
#ls

#mergelist=""
#for f in *${outvcf}
#do
#mergelist+="-V $f "
#done

#echo $mergelist

#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out ${outvcf}

#variants to table
#java $tmj -XX:ParallelGCThreads=1 -Xmx100g -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable \
#-R $ref  \
#-V $outvcf -AMD \
#-o vartable.${outvcf/.vcf.gz/.tsv} \
#-F CHROM -F POS -F QUAL -F QD -F AC -F AN -F DP -F MQ -F MQRankSum -F ReadPosRankSum -F FS -F HET -F SOR

#rm -r `dirname $ref`
#rm  *g.vcf*

cp $SCRATCHDIR/* ${wdir}/${outdir} || export CLEAN_SCRATCH=false
