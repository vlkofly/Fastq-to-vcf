echo  "sample runtime/ncpu readsout totalreads percout duplic percdupli badmate percbadmate badcigar percbadcigar hcm perchcm"
for f in `ls V1*.o*` #change accordingly
do
ft=`tail -n 100 $f`
#echo  $ft
runtime=`grep "ProgressMeter - Total runtime" $f | cut -d"," -f4 | cut -d' ' -f2` # runtime in minutes 
badcigar=`grep "failing BadCigarFilter" $f | cut -d">" -f 2| cut -d' ' -f 2,4|tr -d "("`
dupli=`grep "failing DuplicateReadFilter" $f | cut -d">" -f 2| cut -d' ' -f 2,4|tr -d "("`
badmate=`grep "failing BadMateFilter" $f | cut -d">" -f 2| cut -d' ' -f 2,4|tr -d "("`
hcm=`grep "failing HCMappingQualityFilter" $f | cut -d">" -f 2| cut -d' ' -f 2,4|tr -d "("`
totreads=`grep "reads were filtered out during the traversal out of approximately" $f | cut -d"-" -f 2 | cut -d' ' -f 2,13,16| tr -d "()"`
f=`echo $f| sed 's/\./\t/g' `

echo $f $runtime $totreads $dupli $badmate $badcigar $hcm | tr ' ' \\t
done
