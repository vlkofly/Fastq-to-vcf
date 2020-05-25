# Fastq_to_vcf

**Pipeline to turn whole-genome resequencing data (fastq) to a variant calling format file (vcf)**

The pipeline consists of filtering raw fastq data, mapping them on a reference and calling and filtering genotypes by GATK best practices.

The pipeline is optimised to run on servers with [PBS professional](https://www.altair.com/pbs-works-documentation/) job scheduling system. It can be directly run on [MetaCentrum](https://metavo.metacentrum.cz/en/about/index.html) - Czech computing infrastructure, but on other servers one will have to modify program names, locations, etc.

Major variables for each step of the pipeline can be defined either with an argument or inscript.

0. 0_indexreference.sh: Prepare the reference genome in fasta format - create all the necessary indices for mapping the reads with BWA.
You will need to supply relative path to the fasta file with argument ref. `qsub -v 'ref=yourreference.fasta' 0_indexreference.sh`

1. 1_mapreads.py:  
 

