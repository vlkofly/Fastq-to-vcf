#! /bin/bash

# first argument is a table from parse_variant_table3.R excess_depth_possitions_count.table.tsv
# second argument number of individuals where a site has to be excessive depth

echo "parsing $1 with threshold value $2 individuals"
echo "total number of excess sites"


t=$1
n=$2
wc -l $t

Rscript ~/scripts/hist_n.R $t # print histogram of individual counts
basename=${t/excess_depth_possitions_count.table.tsv/}
echo $basename
cat $t | awk -v n="$n" '{if ($4>n && $2!="NA") {printf "%s\t%s\n", $2, $3  }}' | tee ${basename}.tmp| tr \\t : > ${basename}.${n}_indiv.GATK.intervals # export GATK intervals of DP excess sites in more than N individuals
head ${basename}.tmp

awk '{printf "%s\t%d\t%d\n" ,$1,($2 - 1),$2}' ${basename}.tmp > ${basename}.${n}_indiv.bed #create bed file   

sort-bed ${basename}.${n}_indiv.bed > ${basename}.${n}_indiv.sorted.bed

bedops --ec -m ${basename}.${n}_indiv.sorted.bed | tee ${basename}.${n}_indiv.sorted.merged.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM "sites masked"}'

###use the bed file to check intersection with other mask:
# bedops -i depthmask.allscaff.10_indiv.sorted.merged.bed ../../excess_depth_mask/Blacklist_ExcSites1.6x_maxRD/allscaff.sorted.merged.bed | awk '{SUM += $3-$2} END {print SUM}'


rm ${basename}.tmp

wc -l ${basename}*bed
