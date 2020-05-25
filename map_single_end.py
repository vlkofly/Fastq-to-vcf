#!/usr/bin/env python3

## Usage and help:
# This script takes raw fastq data, checks their quality by fastqc, trims low quality reads by TRIMMOMATIC,
# maps trimmed reads to reference genome while adding readgroup by BWA, removes duplicates by PICARD,
# and finally produces basic statistics on every process that are parsable by multiqc.
# In terms of previous pipeline it consist of steps P2, M1, M2

# This python script is written for usage in metacentrum computing environment and it generates a bash script
# for every sample and sends it to computational queue by command qsub

### Previous version 1.0 is in brno7-cerit, this version is 1.1.
### Changes from previous version:
# Depth statistics computed immediately, and depth file deleted
# argument trial can have parameter r, which will run the rest of the samples after successfull trial run with first sample   
# Input sample names has changes now it has to have two columns first with the old names, the names of files, second new names according majda renaming 
#modified December 2018 for 600_genomes project

import os, sys, argparse, subprocess, glob, re

#create variables that can be entered as an arguments in command line
parser = argparse.ArgumentParser()
parser.add_argument('-samplenames', type=str, metavar='sn', required=True, help='REQUIRED:text file with names of samples e.g 16-AA208n')
parser.add_argument('-datadir',type=str, metavar='datadir', required=True, help='path to the directory where fastq samples are stored')
parser.add_argument('-wd',type=str, metavar='wd', required=True, help='`pwd`')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
parser.add_argument('-trial', type=str, metavar='trial', default='false', help='If changed to true then only one job will be submitted to Metacentrum to check what is going on')
parser.add_argument('-ref', type=str, metavar='ref', required=True, help='full path to reference, the folder must contain indexes and dictionaries e.g. created by indexreference.sh')
args = parser.parse_args()

# input variables
wt = '20:00:00'
ncpu = 4
mem = 90
scratch = 120
adapters = "TruSeq3-SE.fa" # change to SE adapters 
#refpath = "/storage/brno3-cerit/home/filip_kolar/JIC_reference"
refpath = os.path.split(args.ref)[0]
ref = os.path.split(os.path.split(args.ref)[0])[1]+"/"+os.path.split(args.ref)[1] # name of reference genome with the folder

rext = '.fastq.gz' #read extension fg.gz, fastq.gz, .gz
resdir = '/SE_p2m1m2result'
coding = 'JIC_reference/LyV2_genes_only.bed'
tag = "p2m1m2_"
platform = 'illumina'

# create result directory
if os.path.exists(args.wd+resdir) == False:
    os.mkdir(args.wd+resdir)

#read input sample names
listfile = open(args.samplenames,'r')
samples = listfile.readlines() # read the samplenames
listfile.close()

print ("No of samples supplied in the sample name list "+str(len(samples)))

# collect fastqfiles in datadir
files = []
for dp,dirn,filenames in os.walk(args.datadir):
    for filename in [f for f in filenames if f.endswith(rext)]:
        files.append(os.path.join(dp,filename))
#print(files)
print ("Number of files found "+str(len(files)))
# will run one job per sample no matter how many files each sample has


if args.trial == 't': # trial run with only one job/sample
    samples = samples[0:1]

if args.trial == 'r': # run after trial all samples except the firs one, added in version 1.1
    samples = samples[1:len(samples)] # running only first 20 samples

print(samples)

### loop through samples and new names

for s in samples:
    line = s.split("\t")
    s = line[0] 
    n = line[1]
    adapters = line[2].strip("\n")
    rr = (".*"+s+".*")
    r = re.compile(rr)
    f =list(filter(r.match, files)) # search for the sample in the list of files and create list of files corresponding to each individual
    print("this is the new name:"+n+" and this is the old one: "+s)
    nsamp = len(f) # number of files per individual
    f = sorted(f)

    print(s, nsamp)



    sh = open(tag+n+'.sh','w')
    R1path = f[0]

    lanename = R1path.split("/")[1].split('.')[0] # extract lane name from the path datadir/lanenum.blabla/samplefolder/fastq.gz
    Fsname = R1path.split("/")[-1] #extract the filename from path
    logname = 'logstat/trim_'+s+'.log'
    outbase = 'trim_'+Fsname

    ### initialize bash script
    sh.write('#!/bin/bash\n'+
     '#PBS -N '+tag+n+'\n' +
     '#PBS -l walltime='+wt+'\n' +
     '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch)+'gb\n' +
     '#PBS -j oe\n\n')
    ###clean scratch
    sh.write("trap 'clean_scratch' TERM EXIT\n\n")
    sh.write('if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n')
    
    ###copy raw data to scratch
    sh.write('cp '+args.wd+'/'+R1path+' $SCRATCHDIR\n')
    sh.write('cp -r ' + refpath + ' $SCRATCHDIR/\necho input files copied at `date`\n\n')  # load reference to home
    #sh.write('cp '+datadir+adapters+' $SCRATCHDIR\n')
    sh.write('cd $SCRATCHDIR\n')
    # add modules
    sh.write('module add trimmomatic-0.36\n')
    sh.write('module add bwa-0.7.15\n')  # use BWA
    sh.write('module add samtools-1.4\n')  #
    sh.write('module add picard-2.8.1 \n')
    sh.write('module add fastQC-0.11.5 \n\n')

    jav = 'java -XX:ParallelGCThreads=' + str(ncpu-1) + ' -Xmx' + str(mem) + 'g -jar' #java string

    ####make directories
    sh.write('mkdir outtrim\n')
    sh.write('mkdir logstat\n\n')
    
    ####fastqc
    sh.write('fastqc -o logstat -t '+str(ncpu)+' -f fastq '+Fsname+' \n')    

    #####trimminng
    sh.write('cp /software/trimmomatic/0.36/adapters/TruSeq3-SE.fa $SCRATCHDIR\n') # copy adapters # use those default adapters for SE?
    sh.write(jav+' /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar SE' +
    ' -threads '+str(ncpu-1)+' -trimlog '+logname+' '+Fsname+'  '+outbase+' ' +
    'ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:15  SLIDINGWINDOW:4:15 MINLEN:50\n') # trimmomatic job, HEADCROP:8 and MINLEN:75 added
    sh.write('mv trim_* outtrim/\necho trimming finished at `date`\n\n') # move trimmed files to outtrim

    ####bwa mapping
    noextname = Fsname.replace(rext,"") #remove extension from the filename
    sortname = 'sort_'+noextname # add sort

    # you can actually add the readgroup line in the bwa mem already with the option -R
    sh.write('lanename=`zcat '+Fsname+' | head -n 1 | cut -d: -f1,2,3,4`\n')#readgroup name derived from the reads names assume format Casava 1.8
    readgroup = "-R '@RG\\tID:'\"$lanename\"'\\tLB:'\"$lanename\"'one\\tPL:"+platform+"\\tSM:"+n+"'" # here I could play more with the readgroups LB I should do it more properly according how the libs were prepared but in the case of PCR free lib I think it is allright to

    samtools = 'samtools view -@ '+str(ncpu)+' -bu - | samtools sort - -@ '+str(ncpu)+' -m '+str(int((mem-1/6*mem)/ncpu))+'G -o '+sortname

    sh.write('bwa mem -M '+readgroup+' -t '+str(ncpu)+' '+ref+' outtrim/'+outbase+' '+
    ' | '+samtools+'_se.bam || exit 1  \necho reads mapped at `date`\n\n') # mapping of trimmed reads


    ####indexing
    sh.write('samtools index '+sortname+'_se.bam' '\n')

    sh.write(jav+' $PICARD281 BuildBamIndex I='+sortname+'_se.bam\necho all reads sorted and merged `date`\n\n')

    ###deduplicating
    sh.write(jav+' $PICARD281 MarkDuplicates I='+sortname+'_se.bam O='+noextname+'dedup.bam' +
             ' M=logstat/Dup_metrics'+s+'.log ASSUME_SORTED=true TAGGING_POLICY=All\necho deduplicated `date`\n\n')
    sh.write(jav+' $PICARD281 BuildBamIndex I='+noextname+'dedup.bam\necho all reads sorted and merged `date`\n\n')

    ###samtools final stats
    sh.write('samtools idxstats '+noextname+'dedup.bam > logstat/'+s+'.idxstat.txt \n')
    sh.write('samtools flagstat '+noextname+'dedup.bam > logstat/'+s+'.flagstat.txt \n')
    sh.write('samtools stats '+noextname+'dedup.bam > logstat/'+s+'.summary.txt \n')
    sh.write("samtools depth "+noextname+"dedup.bam | awk '{sum+=$3} END {print sum/NR}'" + 
        " > logstat/depth_"+n+".txt \n")
    sh.write("samtools depth -b "+coding+"  "+noextname+"dedup.bam | awk '{sum+=$3} END {print sum/NR}'" +
        " > logstat/codingdepth_"+n+".txt \necho stats done at `date`\n\n") # the file will contain one number: mean of coverage

    ###rename final files, in the multiqc I will have to rename samples
    sh.write('mv '+noextname+'dedup.bam '+n+'_final.bam\n')
    sh.write('mv '+noextname+'dedup.bai '+n+'_final.bai\n\n')
    ###remove unwanted stuff
    sh.write('rm -r `dirname '+ref+'` \n') #remove reference
    sh.write('rm -r outtrim\n')
    sh.write('rm '+Fsname+'\n') # remove raw data
    sh.write('rm *_se.ba*\n\n')
    ###export results
    sh.write('cp -r * '+args.wd+resdir+'/ || export CLEAN_SCRATCH=false\necho results exported at `date`\n')
    sh.close()

    ###check if bash file should be printed or sent to metacentrum
    if args.print == 'false':
        #send by qsub to metacentrum
        cmd = ('qsub '+tag+n+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
        print('\nSubmitted  shell scripts to the cluster!\nFinished!!\n')
    else:
        file = open(tag+n+'.sh','r')
        data = file.read()
        print("bash script follows:\n\n\n")
        print(data)
        print("end of bash script\n\n\n")
        print(samples)




