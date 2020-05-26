# Fastq-to-vcf


**Pipeline to turn whole-genome resequencing data (fastq) to a variant calling format file (vcf)**

The pipeline consists of filtering raw fastq data, mapping them on a reference and calling and filtering genotypes by GATK best practices.
It has been gradually developed by Levi Yant, Jeff DaCosta, Christian Sailer and Jakub Vlček 

The pipeline is optimised to run on servers with [PBS professional](https://www.altair.com/pbs-works-documentation/) job scheduling system. It can be directly run on [MetaCentrum](https://metavo.metacentrum.cz/en/about/index.html) - Czech computing infrastructure, but on other servers one will have to modify program names, locations, etc.

Major variables for each step of the pipeline can be defined either with an argument or inscript.

## 0_indexreference.sh: It prepares the reference genome in fasta format - creates all the necessary indices for mapping the reads with BWA.
You need to supply relative path to the fasta file with argument ref. 
Example how to run the script: `qsub -v 'ref=yourreference.fasta' ~/Fastq-to-vcf/0_indexreference.sh`

## 1_mapreads.py:  This Python script generates separate shell script per every sample that contain commands to: 
- Analyse quality of fastq reads with Fastqc (ver. 0.11.5)
- Trim bad quality reads with TRIMMOMATIC (ver. 0.36) 
- Map reads to the reference with BWA mem (ver. 0.7.15)
- Mark duplicates with Picard (ver. 2.8.1)

Output is sorted and deduplicated bam file, accompanied with multiple descriptive statistics of mapping that can be analysed by a program Multiqc.
You have to supply samplename file that links name of a fastq file with a sample name and adaptors that will be used for trimming: fastq_name\tnew_name\tadaptors\n
Other required arguemnts are -wd = working directory; -datadir = directory with fastq files.
You can also run the script only for the first sample (-trial t) from the list to check errors and if there are none run the analysis for the rest of the samples (-trial r)
Also you can only print the bash script instead of submitting it by qsub if you set -print t
See further options and guidlines within the script. 
Example how to run the script: `python3 ~/Fastq-to-vcf/1_mapreads.py -samplenames sn.tsv -wd "pwd" -datadir ../fastq_data -trial t -ref fullpath_to_reference_folder`

2_callvars.py: This script takes the bam files from previous step and proceeds with per individual variant calling via the GATK tool 'HaplotypeCaller.
You have to supply a -ploidyfile with one column of samplenames and second column with ploidy level 2/4. 
A seperate shell script is created for each sample and sent via qsub to the METACENTRUM cluster. 
If the output directory does not exist, it is created automatically.
You can also run the script only for the first sample (-trial t) from the list to check errors and if there are none run the analysis for the rest of the samples (-trial r)
Also you can only print the bash script instead of submitting it by qsub if you set -print t
See further options and guidlines within the script.
Example how to run the script: `python3 ~/Fastq-to-vcf/2_callvars.py -ploidyfile ploidy.tsv -o ../HC -workdir "pwd"`

3_genotypeGVCF.sh: This script takes the GVCFs files from previous step and creates the final vcf by joint genotyping.
You have to supply following variables
-samples=file containing list of samples
-outvcf=name of output vcf file
-sourcedirs=list of directories (fullpath) with GVCFs files (input files will be fetched from there) 
-outdir=output directory
-ref=full path to a reference folder
-sc=file with scaffolds/chromosomes to parallelize over. One scaffold per line. This allows to run separate jobs for each scaffold.
-NV=yes if you want to run the joint genotype calling for both invariant and variant sites
See further options and notes within the script.
Example how to run the script: `qsub -v 'samples=v3.samplelist.txt,outvcf=outvcf.vcf.gz,sourcedirs=/home/analysis/HC,outdir=finalvcf,ref=fullpath_to_reference_folder' ~/Fastq-to-vcf/3_genotypeGVCF.sh`

 
