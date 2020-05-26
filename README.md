# Fastq-to-vcf


## Pipeline to turn whole-genome resequencing data (fastq) to a variant calling format file (vcf)

The pipeline consists of filtering raw fastq data, mapping them on a reference and calling and filtering genotypes by GATK best practices.
It has been gradually developed by Levi Yant, Jeff DaCosta, Christian Sailer and Jakub Vlček 

The pipeline is optimised to run on servers with [PBS professional](https://www.altair.com/pbs-works-documentation/) job scheduling system. It can be directly run on [MetaCentrum](https://metavo.metacentrum.cz/en/about/index.html) - Czech computing infrastructure, but on other servers one will have to modify program names, locations, etc.

Major variables for each step of the pipeline can be defined either with an argument or inscript.

### 0_indexreference.sh
#### It prepares the reference genome in fasta format - creates all the necessary indices for mapping the reads with BWA.

You need to supply relative path to the fasta file with argument ref. 
Example how to run the script: `qsub -v 'ref=yourreference.fasta' ~/Fastq-to-vcf/0_indexreference.sh`

### 1_mapreads.py
#### This script generates separate shell script per every sample that contains commands to: 
- Analyse quality of fastq reads with Fastqc (ver. 0.11.5)
- Trim bad quality reads with TRIMMOMATIC (ver. 0.36) 
- Map reads to the reference with BWA mem (ver. 0.7.15)
- Mark duplicates with Picard (ver. 2.8.1)

Output is sorted and deduplicated bam file, accompanied with multiple descriptive statistics of mapping that can be analysed by a program Multiqc.  
You have to supply samplename file that links name of a fastq file with a sample name and adaptors that will be used for trimming (fastq_name\tnew_name\tadaptors\n)  
Other required arguemnts are -wd = working directory; -datadir = directory with fastq files.  
You can also run the script only for the first sample (-trial t) from the list to check errors and if there are none run the analysis for the rest of the samples (-trial r)  
Also you can only print the bash script instead of submitting it by qsub if you set -print t  
See further options and guidlines within the script.  
Example how to run the script: 
`python3 ~/Fastq-to-vcf/1_mapreads.py -samplenames sn.tsv -wd "pwd" -datadir ../fastq_data -trial t -ref fullpath_to_reference_folder`

### 2_callvars.py
#### This script takes the bam files from previous step and proceeds with per individual variant calling via the GATK tool 'HaplotypeCaller'.
- All following steps use GATK (ver. 3.7)
 
You have to supply a -ploidyfile with one column of samplenames and second column with ploidy level. (samplename\t2\n)  
A seperate shell script is created for each sample and sent via qsub to the METACENTRUM cluster.   
If the output directory does not exist, it is created automatically.  
You can also run the script only for the first sample (-trial t) from the list to check errors and if there are none run the analysis for the rest of the samples (-trial r)  
Also you can only print the bash script instead of submitting it by qsub if you set -print t  
See further options and guidlines within the script.  
Example how to run the script:  
`python3 ~/Fastq-to-vcf/2_callvars.py -ploidyfile ploidy.tsv -o ../HC -workdir "pwd"`

### 3_genotypeGVCF.sh
#### This script takes the GVCFs files from previous step and creates the final vcf by joint genotyping.
You have to supply following variables  
-samples=file containing list of samples (substring of the GVCFs file) one sample a line  
-outvcf=name of output vcf file  
-sourcedirs=list of directories (fullpath) with GVCFs files (input files will be fetched from there)  
-outdir=output directory 
-ref=full path to a reference folder  
-sc=file with scaffolds/chromosomes to parallelize over. One scaffold per line. This allows to run separate jobs for each scaffold.  
-NV=yes if you want to run the joint genotype calling for both invariant and variant sites  
See further options and notes within the script.  
Example how to run the script:  
`qsub -v 'samples=v3.samplelist.txt,outvcf=outvcf.vcf.gz,sourcedirs=/home/analysis/HC,outdir=finalvcf,ref=fullpath_to_reference_folder' ~/Fastq-to-vcf/3_genotypeGVCF.sh`

### 4_filter.sh
#### This script filters the raw vcf file from previou step according to GATK best practices
You have to supply following variables
-vcf_var= relative path to vcf file with only variants  
-vcf_all (optional) = relative path to vcf with all sites (source of invariant sites, if you want to get number of callable sites)  
-ref=full path to a reference folder  
-outdir=output directory  
-sc=file with scaffolds/chromosomes to parallelize over. One scaffold per line. This allows to run separate jobs for each scaffold.  
-simple=yes if you want to run simplified version of filtering - see details inscript.  
You also have to specify mask files, if you want to remove some sites from vcf (defined within the script) and fourfold site annotation  
Here is the outline of the procedures encoded within the script and description of the output files:  
##### Script outline:
1. Hardfilter biallelic sites (BI) with GATK best practices 
2. Select BI sites that passed the filters
3. Do extra-filtering for genotype missingness with fourfold BI sites
4. Hardfilter non-variant sites (NO) qual threshold 15
5. Select non-variant sites that passed
6. Merge passed BI and NO 
7. Remove the sites specified in masking files Excess heterozygous and excess depth sites (merged.filtered)
8. Get only fourfold sites 
9. Generate variants to table and count number of sites in each step
##### Output files:
- Names derived from input file names, by suffix alternation.
- Vcf files ends with .gz
- Vcf index ends with .tbi (always keep it together with the vcf you use)
- Log file for each step ends with .log
- Key files ready to use in bold, other files kept for debugging and you-never-know-when-you-need purposes.
- Key strings follows consecutively as the corresponding files are created in the process:  
**VCF FILES** 
- novarsel = invariant sites only
- novarfilt = invariant sites after filtration (annotation added)
- novarpassed = invariant sites that have passed the filters
- bisel = biallelic snps only
- bifilt = biallelic snps after filtration
- bipassed = only biallelic snps that have passed filters
- merged = novarpassed and bipassed merged into one vcf 
- **merged.filtered = final vcf with masked sites removed**
- **fourfold.filtered = same as merged.filtered but restricted to fourfold sites**
- fourfold.dp8 = fourfold biallelic snps filtered for depth 8
- fourfold.dp8nc = fourfold biallelic snps filtered for depth 8, GF that did not passed set to nocall
- **fourfold.dp8nc.m0.5 = fourfold biallelic snps with max fraction of missing individual genotypes 0.5**  
**DESCRIPTIVE FILES**
- numbers_in_vcf_files.ssv = number of sites in vcf files
- prefix vartable. = output of command variantsToTable, generated for selection of vcf files and later this table is taken by rscript **parse_variant_table3.R**  that prints some aspects of the results
- filtering_table.pdf/.tsv = what filters caused failing of a variant
- filtering_stats.pdf = distribution of values of filtering parameters, when possible, distribution of passed vs all sites is given
- scatter_depth_chrom.png = depth scatter plot over one scaffold
- boxplot_sampledepth.png = boxplots of sample depths calculated over one scaffold
- boxplot_sampledepthmax100.png = the same as previous with ylim 100 to suppress outliers.
- boxplot_samplequalitymax100.png
- boxplot_samplequality.png
- sample.depth.summary.tsv = depth summary for each sample over all scaffolds
- sample.quality.summary.tsv = genotype quality (GQ) summary for each sample over all scaffolds
- excess_depth_positions_count.table.tsv = positions that showed excessive depth mean+2sd, column n shows in how many samples the limit was exceeded.

Example how to run the script:  
`qsub -v 'vcf_var=snp/halleri_snp_raw.vcf.gz,vcf_all=all/halleri_all_raw.vcf.gz,ref=fullpath_to_reference_folder'  ~/Fastq-to-vcf/4_filter.sh`

### other scripts

#### map_single_end.py
version of 1_mapreads.py script that maps fastq reads from single end sequencing
#### map_unsorted.py
version of 1_mapreads.py script that sorts fastq reads before mapping
#### parse_variant_table3.R
R script to analyse filtering statistics.
It takes one argument, input file, table from command variants to table
#### run_variant_stat.sh
separate script that generates variant to table and then analyse it with parse_variant_table3.R
#### V3_NOparallel.sh
version of 3_genotypeGVCF.sh that does not use parallelisation per scaffold 


