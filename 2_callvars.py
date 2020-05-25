#!/usr/bin/env python3
# This script proceeds with variant calling per individual.

# written in perl by Levi Yant, streamlined by Jeff DaCosta,  by Christian Sailer, 2 May 2016
# adjusted for metacentrum by Jakub Vlcek 5 May 2017 
# adjusted 10. October 2017 ver1.1
# 11. November 2019 add reference option, the same as in 1_mapreads.py

import os, sys, argparse, subprocess, glob

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the analysis ready reads mapped to a reference for variant calling via the GATK tool '+
								 'HaplotypeCaller. You have to supply a file with one column of samplenames and second column of ploidy 2/4. A seperate '+
								 'METACENTRUM shell script is created for each input file and sent to the METACENTRUM cluster. If the output '+
								 'directory does not exist, it is created automatically.')

parser.add_argument('-workdir', type=str, metavar='workdir', help='REQUIRED: Full path to directory with input realigned bam files, use `pwd` ')
parser.add_argument('-ploidyfile', type=str, metavar='diploids_ind', required=True, help='REQUIRED: Full path to text file containing the list of samplenames and ploidy (tab delimited). One individual per line.')
parser.add_argument('-minbaseq', type=int, metavar='min_base_qual', default='20', help='Minimum base quality to be used in analysis [25]')
parser.add_argument('-minmapq', type=int, metavar='min_mapping_qual', default='20', help='Minimum mapping quality to be used in analysis [25]')
parser.add_argument('-suffix', type=str, metavar='HC', default='HC', help='Suffix added to end of output filename [HC1]')
parser.add_argument('-o', type=str, metavar='HC_dir', default='../HC/', help='Realtive path to output directory [HC/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='4', help='Number of requested cores for job [4]')
parser.add_argument('-mem', type=str, metavar='memory', default='16', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='40:00:00', help='Time for job [6-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
parser.add_argument('-trial', type=str, metavar='trial', default='false', help='If changed to true then only one job will be submitted to Metacentrum to check what is going on')
parser.add_argument('-ref', type=str, metavar='ref', required=True, help='full path to reference, the folder must contain indexes and dictionaries e.g. created by indexreference.sh')


args = parser.parse_args()

count = 0
#inscript variables please change according your data
refpath = os.path.split(args.ref)[0]
ref = os.path.split(os.path.split(args.ref)[0])[1]+"/"+os.path.split(args.ref)[1] # name of reference genome with the folder


tag = "final" # a key string to recognize the bam file from other bam files in wd
with open(args.ploidyfile,'r') as p: # in this file there are sample names in firs column and ploidy as integer 2 or 4 in second column sep=\t
	ploidy = p.readlines()

print (str(len(ploidy))+" samples found")
print (ploidy)

#check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
	os.mkdir(args.o)


if args.trial == 't': # trial define shorter list
	ploidy = ploidy[0:1]

if args.trial == 'r': # trial define shorter list
	ploidy = ploidy[1:len(ploidy)]

ncpu=str(args.c)
mem=str(args.mem)
time=str(args.time)
jav = 'java -XX:ParallelGCThreads='+ncpu+' -Xmx'+mem+'g'

#loop through individuals and create one job per individual select ploidy from the ploidy file
count = 0
for s in ploidy:
	sname = s.split("\t")[0]
	sploid = s.split("\t")[1].strip("\n")
	print('\n\n processing sample '+sname+' of ploidy: '+sploid+'\n\n')
	bampath = args.workdir+'/'+sname+'*'+tag+'.ba*'
	bamfile = glob.glob(sname+'*'+tag+'.bam')
	if not bamfile : #added to check if the bamfile is present
	    print ('bam file absent '+sname+'\n')
	    continue
	bamfile = bamfile[0]
	baifile = bamfile.replace('bam','bai')
	print (bamfile)
	print (baifile)
	shf = open(args.o+'V1'+sname+'.sh','w')
	#write METACENTRUM shell file
	shf.write('#!/bin/bash \n'+
			'#PBS -N V1.'+sname+'\n'+
			'#PBS -l walltime='+time+'\n' +
			'#PBS -l select=1:ncpus='+ncpu+':mem='+mem+'gb:scratch_local=40gb\n' + 
			'#PBS -j oe\n\n'
			'module add gatk-3.7 \n'+
			'module add samtools-1.4 \n'+
			'trap \'clean_scratch\' TERM EXIT\n'+
			'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
			'cp -r '+refpath+' $SCRATCHDIR\n' +
			'cp '+bampath+' $SCRATCHDIR\n'+
			'cd $SCRATCHDIR\necho data loaded at `date`\n\n')
	
	###index bam file if index got lost
	#shf.write('samtools index '+bamfile+' '+baifile+' \n' )
	#HaplotypeCaller string
	shf.write(jav+' -jar $GATK/GenomeAnalysisTK.jar -T '+
		'HaplotypeCaller -I '+bamfile+' --min_base_quality_score '+str(args.minbaseq)+' --min_mapping_quality_score '+
		str(args.minmapq)+' -rf BadMate '+
		'-R '+ref+' -o '+sname+'_'+str(args.suffix)+'.g.vcf.gz -ploidy '+sploid+' -stand_call_conf 20 -ERC GVCF '+
		'--pcr_indel_model NONE -nct '+ncpu+' --max_num_PL_values 350 \necho calling done at `date`\n\n')
			
	
	shf.write('rm -r '+ref+'\n'+
		'rm '+sname+'*'+tag+'.ba* \n'+ 
		'cp *vcf* '+args.workdir+'/'+args.o+' || export CLEAN_SCRATCH=false\n'+
		'printf "\\nFinished\\n\\n"\n')
	shf.close()

	#check if METACENTRUM shell file should be printed or sent to METACENTRUM 
	if args.print == 'false':
		#send METACENTRUM job to METACENTRUM cluster
		cmd = ('qsub '+args.o+'V1'+sname+'.sh')
		p = subprocess.Popen(cmd, shell=True)
		sts = os.waitpid(p.pid, 0)[1]
	else:
		file = open(args.o+'V1'+sname+'.sh','r')
		data = file.read()
		print(data)

	count +=1
#if appropriate, report how many METACENTRUM shell files were sent to METACENTRUM 
if args.print == 'false':
	print('\nSent '+str(count)+' jobs to the METACENTRUM\n\nFinished!!\n\n')
