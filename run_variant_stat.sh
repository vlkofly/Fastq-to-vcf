#!/bin/bash
#PBS -N parse_vartable 
#PBS -l walltime=60:00:00
#PBS -l select=1:ncpus=1:mem=350gb:scratch_local=300gb
#PBS -j oe

# run parse_variant_table for tab

if [ -z "$tab" ]; then
	echo "Error! supply result of variantToTable"
	echo "qsub -v 'tab=blabla.tsv'"
	exit 1
fi

wd=$PBS_O_WORKDIR
script="/storage/brno3-cerit/home/vlkofly/scripts/parse_variant_table3.R"
base=${tab/.tsv/}
nov=n # for variants 
#nov=novar # uncomment for non-variants
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp $wd/$tab $SCRATCHDIR/

cd $SCRATCHDIR

Rscript $script $tab $nov

rm $tab

mkdir $base
mv * $base/
cp -r $base $wd/ || export CLEAN_SCRATCH=false

