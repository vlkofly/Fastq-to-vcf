#!/bin/bash
#PBS -N F1_filter
#PBS -l walltime=256:00:00
#PBS -l select=1:ncpus=8:mem=650gb:scratch_local=800gb
#PBS -j oe

# This script takes the results of genotypeGVCF V3 filters data according GATK best practices.
# It can filter either only varible sites or both variable and invariable sites (to do both one must supply vcf_var and vcf_all - vcf with only variants an vcf with all sites)
# To run it with your data you need to supply following informations:


if [ -z "$vcf_var" ]; then
	echo "Error! Specify relative path (relative to working directory) to vcf file after joinGVCF calling that contains only variants!"
        echo "qsub -v 'vcf_var=blablabla.vcf.gz'"
        exit 1
        fi

if [ -z "$vcf_all" ]; then
        echo "Vcf with all sites not supplied, only snps will be filtered!"
	vcf_all=no
        echo "qsub -v 'vcf_all=blablabla.vcf.gz'"
        fi
if [ -z "$outdir" ]; then
        echo "default outdir used"
	outdir=filtering_results
        fi

if [ -z "$ref" ]; then
        echo "Specify reference full path, default reference used instead"
        echo "qsub -v '/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta'"
        ref=/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta
        fi

if [ -z "$sc" ]; then
        echo "Specify full path to scaffold list to run across in parallel"
        echo "if not specified default for lyrata is used"
        echo "one scaffold/chromosome per line"
        sc=/storage/brno3-cerit/home/filip_kolar/lyrata_parallel_scaffold_list.txt
        fi

if [  "$simple" == "yes" ]; then
        echo "simplified version of filtering will be run"
        echo "without fourfold sites filtering and masking"
	echo "used e.g. for plastome"
else
	echo "full version filtering run"
	fi

echo "started at `date`"

####define some other variables ( those variables won't change much so I define them inscript)

wd=$PBS_O_WORKDIR
nt=8 # number of threads change based on number of scaffolds - but it is better to run one job per a scaffold

refdir=`dirname $ref`
ref=`echo $ref | sed 's/^.*\/\(.*\/.*\)$/\1/g'`
scaf=`basename $sc`

# annotations for masking and fourfold sites
#hetmask="/storage/brno3-cerit/home/filip_kolar/lyrata_annotations/Blacklist.Excess.Het.List.Genes.With.5orMore.2PopsMin.AllHet.Sites.GATK.intervals"
hetmask="/storage/brno3-cerit/home/filip_kolar/halleri_fixedHetPerGene_2orMorePop_5orMoreSites.intervals"
#dpmask="/storage/brno3-cerit/home/filip_kolar/lyrata_annotations/Blacklist_allscaff_ExcSites1.6x_maxRD.GATK.intervals"
dpmask="/storage/brno3-cerit/home/filip_kolar/depthmask.negative.halleri.allscaff.10_indiv.GATK.intervals"
fourfold="/storage/brno3-cerit/home/filip_kolar/lyrata_annotations/four.fold.degenerate.sites.GATK.intervals"


trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

if [ ! -d "${wd}/${outdir}" ] ; then echo "output dir created";  mkdir ${wd}/${outdir}; fi # create output directory if not there yet

###copy everything necessary to SCRATCHDIR
cp -r $refdir $SCRATCHDIR/
cp $wd/$vcf_all $SCRATCHDIR/
cp $wd/${vcf_all}.tbi $SCRATCHDIR/ # copy index as well
cp $wd/$vcf_var $SCRATCHDIR/
cp $wd/${vcf_var}.tbi $SCRATCHDIR/
cp $hetmask $SCRATCHDIR/
cp $dpmask $SCRATCHDIR/ # solve mask from F2
cp $fourfold $SCRATCHDIR/
cp $sc $SCRATCHDIR/ # list of scaffolds to parallelize over
#cp /storage/brno3-cerit/home/filip_kolar/lyrata_scaffolds_rest.intervals $SCRATCHDIR/ # batch of short scaffolds, change with reference TODO 

hetmask=`basename $hetmask` # after copying the data remove the path from annotation variables
dpmask=`basename $dpmask`
fourfold=`basename $fourfold`
vcf_var=`basename $vcf_var`
vcf_all=`basename $vcf_all`

###define output and auxiliary file names
bisnp_sel=${vcf_var/vcf.gz/bisel.vcf.gz} # selected biallelic snps
bisnp_filt=${vcf_var/vcf.gz/bifilt.vcf.gz} # filtered biallelic snps
bisnp_passed=${vcf_var/vcf.gz/bipassed.vcf.gz} # sites that have not passed removed
#novar
novar_sel=${vcf_all/vcf.gz/novarsel.vcf.gz} # selected novariants
novar_filt=${vcf_all/vcf.gz/novarfilt.vcf.gz} # filtered novariants
novar_passed=${vcf_all/vcf.gz/novarpass.vcf.gz} # sites that have not passed filters removed
#merged
merged_novar_bi_passed=${vcf_var/vcf.gz/merged.vcf.gz} # merged biallelic snps and no variants
merged_fourfold=${vcf_var/vcf.gz/snps.fourfold.vcf.gz} # extracted 4 fold degenerated snps after filtering
merged_filtered=${vcf_var/vcf.gz/merged.filtered.vcf.gz} # merged files after masking with 300 masks
merged_fourfold_filtered=${vcf_var/vcf.gz/fourfold.filtered.vcf.gz} # extracted 4 fold degenerated sites, after merging filtering
cd $SCRATCHDIR


#head -n 8 $scaf > lyrata_parallel_scaffold_list_eight.txt # select only subset of scaffolds

echo "data imported succesfully at `date`:"

ls

#load modules (programs)
module add gatk-3.7
module add bcftools-1.6
module add R-3.4.3-gcc
module add picard-2.8.1
module add parallel-20160622


tmj="-Djava.io.tmpdir=$SCRATCHDIR"




###export only biallelic snp sites######
# the scatter gather approach in each of the pipeline step, due to speed improvement,
# it means that command is run separately for each chunk of a genome,
# and output for separated chunks is then merged together

parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants -R $ref -V $vcf_var -L {} -o {.}$bisnp_sel --restrictAllelesTo BIALLELIC -selectType SNP -nt 1 &> {.}${bisnp_sel/.vcf.gz/}.log"

# if you want to merge the intermediate vcfs
#echo "biallelic snps selected in parallel at `date`"
#mergelist=""
#for f in *bisel.vcf.gz
#do
#mergelist+="-V $f "
#done
#
#echo $mergelist
#
#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $bisnp_sel
#
#parallel -a $scaf -j $((nt-1)) "rm {.}$bisnp_sel" # intermedate scattered files removed 
#parallel -a $scaf -j $((nt-1)) "rm {.}${bisnp_sel}.tbi"
#
#
#echo "biallelic snps concatenated at `date` "
#ls

###hard filter the biallelic snps with best practices values

parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $ref \
 -V {.}$bisnp_sel -o {.}$bisnp_filt \
-filter 'QD < 2.0' --filterName 'QD' -filter 'FS > 60.0' --filterName 'FS' \
-filter 'MQ < 40.0' --filterName 'MQ' -filter 'MQRankSum < -12.5' --filterName 'MQRS' \
-filter 'ReadPosRankSum < -8.0' --filterName 'RPRS' -filter 'SOR > 3.0' --filterName 'SOR' -nt 1 &> {.}${bisnp_sel/gz/filtering.log}"

echo "bisnp filtered in parallel `date`"
#mergelist=""

#for f in *bifilt.vcf.gz
#do
#mergelist+="-V $f "
#done

#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $bisnp_filt
#
#
#parallel -a $scaf -j $((nt-1)) "rm {.}$bisnp_filt"
#parallel -a $scaf -j $((nt-1)) "rm {.}${bisnp_filt}.tbi"
#
#echo "biallelic snps filtered and concatenated  at `date`"
#ls

###select only biallelic variants that passed the filters
parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $ref \
      	-V {.}$bisnp_filt  -o {.}$bisnp_passed --excludeFiltered  -nt 1 &> {.}${bisnp_passed}.log"

echo "biallelic snps passing filter exported at `date`"
#ls

##select four fold degenerated snps from bifilt and filter - this is nut run if simple == true
if [ "$simple" != "yes" ]; then
	echo starting with fourfold selection
	parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar \
	 -T SelectVariants -R $ref \
	 -V {.}$bisnp_passed -o {.}$merged_fourfold -nt 1 -XL $hetmask --interval_set_rule INTERSECTION \
	 -XL $dpmask -L $fourfold &> {.}${merged_fourfold}.log"
	
	#rm $merged_fourfold
	
	#mergelist=""
	#for f in *snps.fourfold.vcf.gz
	#do
	#mergelist+="-V $f "
	#done
	
	#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
	#org.broadinstitute.gatk.tools.CatVariants -R $ref \
	#$mergelist -assumeSorted -out $merged_fourfold
	
	#parallel -a $scaf -j $((nt-1)) "rm {.}$merged_fourfold"
	#parallel -a $scaf -j $((nt-1)) "rm {.}${merged_fourfold}.tbi"
	
	echo bisnp fourfold exported in parallel and merged at `date`
	
	### depth and missingness filtering of fourfold degenerated snps
	
	parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V {.}$merged_fourfold \
	--genotypeFilterExpression 'DP < 8' --genotypeFilterName 'DP' -log {.}${merged_fourfold/fourfold.filtered.vcf.gz/fourfold.dp8.vcf.gz}.log \
	 -o {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8.vcf.gz} &>/dev/null"
	
	parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R \
	$ref -V {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8.vcf.gz} -log {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8nc.vcf.gz}.log \
	--setFilteredGtToNocall -o {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8nc.vcf.gz} &> /dev/null" 
	
	
	parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R \
	$ref -V {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8nc.vcf.gz} --maxNOCALLfraction 0.5 \
	-log {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8nc.m0.5.vcf.gz}.log \
	 -o {.}${merged_fourfold/fourfold.vcf.gz/fourfold.dp8nc.m0.5.vcf.gz} &> /dev/null" # the only file from this procedure will be the selected variants with less than 50% missingness
	
	#rm ${merged_fourfold/fourfold.gz/fourfold.dp8.gz}*
	#rm ${merged_fourfold/fourfold.gz/fourfold.dp8nc.gz}*
fi
# if there is no variant file then stop here and export snp results:

if [ "$vcf_all" == "no" ]; then
echo "all sites not supplied, only snps were filtered!"
# prepare depth filtered vcf 
bipas=`ls *bipassed.vcf.gz`
java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V $bipas \
        --genotypeFilterExpression "DP < 8" --genotypeFilterName "DP" \
         -o ${bipas/bipassed.vcf.gz/bipassed.dp8.vcf.gz}

java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R \
        $ref -V ${bipas/bipassed.vcf.gz/bipassed.dp8.vcf.gz}  \
        --setFilteredGtToNocall -o ${bipas/bipassed.vcf.gz/bipassed.dp8nc.vcf.gz}


java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R \
        $ref -V ${bipas/bipassed.vcf.gz/bipassed.dp8nc.vcf.gz} --maxNOCALLfraction 0.5 \
        -o ${bipas/bipassed.vcf.gz/bipassed.dp8nc.m0.5.vcf.gz}

rm ${bipas/bipassed.vcf.gz/bipassed.dp8.vcf.gz}*
rm ${bipas/bipassed.vcf.gz/bipassed.dp8nc.vcf.gz}*

ls *vcf.gz | parallel -j $((nt-1)) 'n=`zcat {} | grep -v "#" | wc -l`; echo {} $n' >> ${vcf_var/vcf.gz/numbers_in_vcf_files.ssv}
parallel -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -Xmx300g -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable \
	-R $ref --showFiltered -log {}vartable.log \
	-V {} -AMD \
	-o vartable.{.}.tsv \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F ReadPosRankSum -F FS -F HET -F SOR -F FILTER \
	-GF DP -GF GQ &>/dev/null" ::: $bisnp_filt
rm ${vcf_var}
rm ${vcf_var}.tbi
find -maxdepth 1 -type f -exec cp {} ${wd}/${outdir}/ \; || export CLEAN_SCRATCH=false
exit 1
fi

#####export only non_variant sites#####
parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $ref \
 -V $vcf_all -L {} -o {.}$novar_sel -selectType NO_VARIATION   -nt 1 &> {.}${novar_sel}.log"

echo "non-variants exported in parallel `date`"
#mergelist=""
#for f in *novarsel.vcf.gz
#do
#mergelist+="-V $f "
#done
#
#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $novar_sel

#parallel -a $scaf -j $((nt-1)) "rm {.}$novar_sel"
#parallel -a $scaf -j $((nt-1)) "rm {.}${novar_sel}.tbi"

#echo "non-variant sites concatenated and at `date`"
ls
###hard filter non variant sites
parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $ref \
 -V {.}$novar_sel -o {.}$novar_filt  \
-filter 'QUAL < 15.0' --filterName 'QUAL'  -nt 1 &> {.}${novar_filt}.log" # QUAL in non-variant is a phred of prob of the site being variant

echo "non-variant sites filtered in parallel at `date`"

#mergelist=""
#for f in *novarfilt.vcf.gz
#do
#mergelist+="-V $f "
#done

#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $novar_filt

#parallel -a $scaf -j $((nt-1)) "rm {.}$novar_filt"
#parallel -a $scaf -j $((nt-1)) "rm {.}${novar_filt}.tbi"

#echo "non-variant sites filtered merged at `date`"
ls
###select only non_variants that passed the filters why this is not parallelised?
parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $ref \
-V {.}$novar_filt  -o {.}$novar_passed --excludeFiltered -nt 1 &> {.}${novar_passed}.log"

echo "non-variant sites passing filter exported at `date`"
ls
#####merge non-variant and biallelic snps

parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar \
  -T CombineVariants \
  -R $ref  \
  --variant:novar {.}$novar_passed \
  --variant:snps {.}$bisnp_passed \
  -o {.}$merged_novar_bi_passed \
  -genotypeMergeOptions PRIORITIZE \
  -priority novar,snps &> {.}${merged_novar_bi_passed}.log"

echo "non-variant sites and biallelic snps merged in parallel at `date`"

#mergelist=""
#for f in *merged.vcf.gz
#do
#mergelist+="-V $f "
#done

#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $merged_novar_bi_passed

#keep the separate chroms
#parallel -a $scaf -j $((nt-1)) "rm {.}$merged_novar_bi_passed"
#parallel -a $scaf -j $((nt-1)) "rm {.}${merged_novar_bi_passed}.tbi"

echo "non-variants and biallelic snps merged and parallel files merged at `date`"
ls

##here only the eight scaffolds are being analysed so the parallel file must be modified

#head -n 8 $scaf > lyrata_parallel_scaffold_list_eight.txt 
#
####masking variants
parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $ref \
 -V {.}$merged_novar_bi_passed -nt 1 -L {} -o {.}$merged_filtered \
 -XL $hetmask -XL $dpmask --interval_set_rule INTERSECTION \
  &> {.}${merged_filtered}.log" # create your own depth mask? The interval file with excess depth sites is only for 8 scaffolds

echo "heterozygosity and depth masked in parallel at `date`"

#rm $merged_filtered

#mergelist=""
#for f in *merged.filtered.vcf.gz
#do
#mergelist+="-V $f "
#done

#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $merged_filtered

#parallel -a lyrata_parallel_scaffold_list_eight.txt -j $((nt-1)) "rm {.}$merged_filtered"
#parallel -a lyrata_parallel_scaffold_list_eight.txt -j $((nt-1)) "rm {.}${merged_filtered}.tbi"

echo "heterozygosity and depth masked and merged at `date`" 
#
ls
####continue with the paralelization here


###select fourfold degen sites from merged_filtered
parallel -a $scaf -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -jar $GATK/GenomeAnalysisTK.jar \
 -T SelectVariants -R $ref \
 -V {.}$merged_filtered -L {} -o {.}$merged_fourfold_filtered \
-nt 1 -L $fourfold --interval_set_rule INTERSECTION -log {.}${merged_fourfold_filtered}.log &>/dev/null"

#rm $merged_fourfold_filtered

#mergelist=""
#for f in *fourfold.filtered.vcf.gz
#do
#mergelist+="-V $f "
#done


#java $tmj -XX:ParallelGCThreads=1 -cp $GATK/GenomeAnalysisTK.jar \
#org.broadinstitute.gatk.tools.CatVariants -R $ref \
#$mergelist -assumeSorted -out $merged_fourfold_filtered

#parallel -a lyrata_parallel_scaffold_list_eight.txt -j $((nt-1)) "rm {.}$merged_fourfold_filtered"
#parallel -a lyrata_parallel_scaffold_list_eight.txt -j $((nt-1)) "rm {.}${merged_fourfold_filtered}.tbi"

echo "fourfold degenerated sites exported at `date`"



####stats
echo "vcf_file numvariants"
ls *vcf.gz | parallel -j $((nt-1)) 'n=`zcat {} | grep -v "#" | wc -l`; echo {} $n' >> ${vcf_var/vcf.gz/numbers_in_vcf_files.ssv}


echo "vartotable processing"
parallel -j $((nt-1)) "java $tmj -XX:ParallelGCThreads=1 -Xmx300g -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable \
-R $ref --showFiltered -log {}vartable.log \
-V {} -AMD \
-o vartable.{.}.tsv \
-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F ReadPosRankSum -F FS -F HET -F SOR -F FILTER \
-GF DP -GF GQ &>/dev/null" ::: scaffold_1lyrata_snp_raw.bifilt.vcf.gz # also fix

for v in scaffold_1lyrata_all_raw.novarfilt.vcf.gz # fix 
do
java $tmj -XX:ParallelGCThreads=1 -Xmx300g -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable \
-R $ref --showFiltered -log ${v}vartable.log \
-V $v -AMD \
-o vartable.${v/.gz/.tsv} \
-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FILTER -F MQRankSum -F ReadPosRankSum -F FS -F HET -F SOR \
-GF DP -GF RGQ &>/dev/null
done

echo "variant tables generated at `date`"
ls *vartable

####bcftools stats
bcftools stats scaffold_1lyrata_snp_raw.fourfold.filtered.vcf.gz  >  scaffold_1lyrata_snp_raw.fourfold.filtered.vcf.bcfstat.txt

#vcftools stats, depth and so on
####stats with R
#cp -r * ${wdf}${outdir}/


# the rstats have to be done separately
#echo 'processing variant table with R'
#again do it in loop for fourfold, bifiltr, nonfiltr, merged_nobi passed

#mkdir R_results
#cd R_results
#cp ${wdf}/scripts/parse_variant_table3.R .
#for v in ../vartable*
#do
#Rscript parse_variant_table3.R $v
#done

#cd ../
#echo "rstats finished at `date`"
#ls R_results

rm $vcf
rm -r JIC_reference # todo universalize 
rm $hetmask
rm $fourfold 
rm $dpmask
cp * ${wd}/${outdir}/ || export CLEAN_SCRATCH=false

echo "filtering of done and results exported at `date`"

