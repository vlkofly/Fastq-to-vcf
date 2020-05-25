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
# >version mocker1
####another version where I will add two more lanes per sample may 2018


#####in October 2018 script p2m1m2_melanotis.py modified to be used for 600_genomes project

### implement renaming again
### more efficient data collection, rather than supplying table with specification of datadir I will put all raw fastq files into one directory.

import os, sys, argparse, subprocess, glob, re

from collections import defaultdict
#create variables that can be entered as an arguments in command line
parser = argparse.ArgumentParser()
parser.add_argument('-samplenames', type=str, metavar='sn', required=True, help='REQUIRED:text file defining fastqname\tnewname\tadapter')
parser.add_argument('-datadir',type=str, metavar='datadir', required=True, help='path to the directory where fastq samples are stored')
parser.add_argument('-wd',type=str, metavar='wd', required=True, help='`pwd`')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
parser.add_argument('-trial', type=str, metavar='trial', default='false', help='If changed to true then only one job will be submitted to Metacentrum to check what is going on')
args = parser.parse_args()

# input variables
wt = '30:00:00'
ncpu = 6
mem = 180 
scratch = 220
adapters = "TruSeq3-PE-2.fa" # from TRIMMOMATIC default, 
refpath = "/storage/brno3-cerit/home/filip_kolar/JIC_reference"
rext = '.fastq.gz' # check for multiple fastq extensions

ref = 'JIC_reference/alygenomes.fasta' # name of reference genome 
platform = 'illumina'
resdir = '/serp_p2m1m2result'
coding = 'JIC_reference/LyV2_genes_only.bed'
tag = 'serp_p2m1m2'

# create result directory
if os.path.exists(args.wd+resdir) == False:
    os.mkdir(args.wd+resdir)

#read input sample names
listfile = open(args.samplenames,'r')
samples = listfile.readlines() # read the samplenames
listfile.close()

print ("No of samples supplied in the sample name list "+str(len(samples)))

# collect all fastqfiles in the datadir
files = []
for dp,dirn,filenames in os.walk(args.datadir):
    for filename in [f for f in filenames if f.endswith(rext)]:
        files.append(os.path.join(dp,filename))

#files = glob.glob(args.datadir+"/*/*"+rext)# this assumes that fastq files are datadir/lanefolder/samplefolder/fastq.gz # change it for polyglottos

print(files)
print ("Number of files found "+str(len(files)))
# will run one job per sample no matter how many files each sample has


if args.trial == 't': # trial run with only one job/sample
     samples = samples[0:1]

if args.trial == 'r': # run after trial all samples except the firs one, added in version 1.1
     samples = samples[1:len(samples)]



### loop through samples and new names ###the new name function deprecated and added again add one more column with adapter information

for s in samples:
    line = s.split("\t")
    s = line[0]
    nw = line[1]
    adapters = line[2].strip("\n")
    rr = (".*"+s+".*") # like this it would even recognize sample if the name is only in the folder
    r = re.compile(rr)
    f =list(filter(r.match, files)) # search for the sample in the list of files and create list of files corresponding to each individual
    nsamp = len(f) # number of files per individual
    print("this is the new name:"+nw+" and this is the old one: "+s+"\nfastq files:\n")
    f = sorted(f)
    print (f)
    print(s, nsamp)
    nlane = int(nsamp/2)
    print ("number of libraries in this sample "+str(nlane)+"\n")


    sh = open(tag+nw+'.sh','w')
    #create a dictionary with lane as a key and its reads as two values
    #lanes = {n: f[n:n+2] for n in range(nlane)} 
    x=0
    lanes = defaultdict(list)
    for n in range(nlane):
        r = f[x:x+2] # it assumes that lanes go after each other
        x+=2
        lanes[n].append(r)
    print(lanes)

    ### initialize bash script
    sh.write('#!/bin/bash\n'+
     '#PBS -N  '+tag+nw+'\n' +
     '#PBS -l walltime='+wt+'\n' +
     '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch)+'gb\n' +
     '#PBS -j oe\n\n')
    ###clean scratch
    sh.write("trap 'clean_scratch' TERM EXIT\n\n")
    sh.write('if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n')
    
    ###copy raw data to scratch
    for l in lanes:
        print(lanes[l][0][0])
        sh.write('cp '+args.wd+'/'+lanes[l][0][0]+' $SCRATCHDIR\n')
        sh.write('cp '+args.wd+'/'+lanes[l][0][1]+' $SCRATCHDIR\n')
    sh.write('cp -r ' + refpath + ' $SCRATCHDIR/\necho input files copied at `date`\n\n')  # load reference to home
    #sh.write('cp '+args.wd+resdir+'/'+s+'merged.ba* $SCRATCHDIR\n')
    sh.write('cd $SCRATCHDIR\n')
    # add modules
    sh.write('module add trimmomatic-0.36\n')
    sh.write('module add bwa-0.7.15\n')  # use BWA
    sh.write('module add samtools-1.4\n')  #
    sh.write('module add picard-2.8.1 \n')
    sh.write('module add parallel-20160622\n')
    sh.write('module add fastQC-0.11.5 \n')
    sh.write('module add parallel-20160622\n\n') 
           
            

    jav = 'java -XX:ParallelGCThreads=' + str(ncpu) + ' -Xmx' + str(mem) + 'g -jar' #java string

    ####make directories
    sh.write('mkdir outtrim\n')
    sh.write('mkdir logstat\n\n')

    ####fastqc of the raw data

    sh.write('ls *'+rext+' | parallel -j '+str(int(ncpu)-1)+' "fastqc -o logstat -t 2 -f fastq {}" \n\n')    
    #####trimminng
    
    sh.write('cp /storage/brno3-cerit/home/filip_kolar/600_genomes/adapters/'+adapters+' $SCRATCHDIR\n\n') # copy adapters
    sh.write('cp /storage/brno3-cerit/home/filip_kolar/fastqCombinePairedEnd.py $SCRATCHDIR\n\n') # copy script that sorts reads
    for l in lanes:
        Fsname = lanes[l][0][0].split("/")[-1] #extract the filename from path
        Rsname = lanes[l][0][1].split("/")[-1]
        logname = 'logstat/trim_lane'+str(l)+s+'.log'
        outbase = 'trim_'+Fsname
        sh.write('python fastqCombinePairedEnd.py '+Fsname+' '+Rsname+' None \n')
        sh.write('mv '+Fsname+'_pairs_R1.fastq.gz '+Fsname+' \n')
        sh.write('mv '+Rsname+'_pairs_R2.fastq.gz '+Rsname+' \n') #owerwrite former fastqfiles 

        #I could also do it with parallel, which would be definitely more efficient. this is the next step
        sh.write(jav+' /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar PE' +
        ' -threads '+str(ncpu-1)+' -trimlog '+logname+' '+Fsname+' '+Rsname+' -baseout '+outbase+' ' +
        'ILLUMINACLIP:'+adapters+':2:30:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50\n\n') # trimmomatic job
        sh.write('mv trim_* outtrim/\necho trimming of lane '+str(l)+' finished at `date`\n\n') # move trimmed files to outtrim
    ####fastqc of trimmed data
    sh.write('ls outtrim/*'+rext+' | parallel -j '+str(int(ncpu)-1)+' "fastqc -o logstat -t 2 -f fastq {}" \n\n')

    ####bwa mapping
    for l in lanes:
        Fsname = lanes[l][0][0].split("/")[-1]
        noextname = Fsname.replace(rext,"") #remove extension from the filename
        sortname = 'sort_'+noextname # add sort it assumes that the name of fastq files for each lane is different

        # you can actually add the readgroup line in the bwa mem already with the option -R
        sh.write('echo "processing lane '+str(l)+'"  \nlanename=`zcat '+Fsname+' | head -n 1 | cut -d: -f1,2,3,4`\n')#readgroup name derived from the reads names assume format Casava 1.8
        readgroup = "-R '@RG\\tID:'\"$lanename\"'\\tLB:'\"$lanename\"'one\\tPL:"+platform+"\\tSM:"+nw+"'" # here I could play more with the readgroups LB I should do it more properly according how the libs were prepared but in the case of PCR free lib I think it is allright to

        samtools = 'samtools view -@ '+str(ncpu)+' -bu - | samtools sort - -@ '+str(ncpu)+' -m '+str(int((mem-1/6*mem)/ncpu))+'G -o '+sortname

   
        sh.write('bwa mem -M '+readgroup+' -t '+str(ncpu)+' '+ref+' outtrim/trim_*'+noextname+'*_1P'+rext+' '+
        'outtrim/trim_*'+noextname+'*_2P'+rext+' | '+samtools+'_paired.bam || exit 1  \necho paired reads mapped at `date`\n\n') # mapping paired reads

        sh.write('bwa mem -M '+readgroup+' -t '+str(ncpu)+' '+ref+' outtrim/trim_*'+noextname+'*_1U'+rext+' '+
        '| '+samtools+'_unpaired1.bam || exit 1 \necho unpaired reads1 mapped at `date`\n\n') # mapping unpaired 1 reads

        sh.write('bwa mem -M '+readgroup+' -t '+str(ncpu)+' '+ref+' outtrim/trim_*'+noextname+'*_2U'+rext+' '+
        '| '+samtools+'_unpaired2.bam || exit 1 \necho unpaired reads2 lane $l mapped at `date`\n\n') # mapping unpaired 2 reads

        ####indexing
        sh.write('samtools index '+sortname+'_paired.bam' '\n')
        sh.write('samtools index '+sortname+'_unpaired1.bam' '\n')
        sh.write('samtools index '+sortname+'_unpaired2.bam' '\necho all reads of lane $l sorted `date`\n\n\n')
    

    ####merge paired and unpaired reads and index
    sh.write('for b in *bam\ndo\nblist+="I=$b "\n done \n echo $blist\n\n')
    sh.write(jav+' $PICARD281 MergeSamFiles $blist O='+s+'merged.bam \n') 

    sh.write(jav+' $PICARD281 BuildBamIndex I='+s+'merged.bam\necho all reads sorted and merged `date`\n\n')

    ###deduplicating
    sh.write("limit=$(echo `ulimit -n` - 50 | bc)\n")
    sh.write(jav+' $PICARD281 MarkDuplicates I='+s+'merged.bam O='+s+'dedup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$limit' +
             ' M=logstat/Dup_metrics'+s+'.log ASSUME_SORTED=true TAGGING_POLICY=All\necho deduplicated `date`\n\n')
    sh.write(jav+' $PICARD281 BuildBamIndex I='+s+'dedup.bam\necho all reads sorted and merged `date`\n\n')

    ###samtools final stats
    sh.write('parallel -j '+str(ncpu-1)+' "samtools {} '+s+'dedup.bam > logstat/'+s+'.{}.txt" ::: idxstats flagstat stats\n\n')

    #sh.write('samtools idxstats '+s+'dedup.bam > logstat/'+s+'.idxstat.txt \n')
    #sh.write('samtools flagstat '+s+'dedup.bam > logstat/'+s+'.flagstat.txt \n')
    #sh.write('samtools stats '+s+'dedup.bam > logstat/'+s+'.summary.txt \n')
    sh.write("samtools depth "+s+"dedup.bam | awk '{sum+=$3} END {print sum/NR}'" + 
        " > logstat/depth_"+s+".txt \n")
    #sh.write("samtools depth -b "+coding+"  "+noextname+"dedup.bam | awk '{sum+=$3} END {print sum/NR}'" +
    #   " > logstat/codingdepth_"+s+".txt \necho stats done at `date`\n\n") # the file will contain one number: mean of coverage
    # rename to new name
    sh.write('mv '+s+'dedup.bam '+nw+'_final.bam\n')
    sh.write('mv '+s+'dedup.bai '+nw+'_final.bai\n\n')


    ###remove unwanted stuff
    sh.write('rm -r `dirname '+ref+'` \n') #remove reference
    #sh.write('rm -r outtrim\n')
    sh.write('rm *paired*.bam*\n')
    sh.write('rm *'+rext+'\n') # remove raw data
    #sh.write('rm *merged.ba*\n\n')
    ###export results
    sh.write('cp -r * '+args.wd+resdir+'/ || export CLEAN_SCRATCH=false\necho results exported at `date`\n')
    sh.close()

    ###check if bash file should be printed or sent to metacentrum
    if args.print == 'false':
        #send by qsub to metacentrum
        cmd = ('qsub '+tag+nw+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
        print('\nSubmitted  shell scripts to the cluster!\nFinished!!\n')
    else:
        file = open(tag+nw+'.sh','r')
        data = file.read()
        print("bash script follows:\n\n\n")
        print(data)
        print("end of bash script\n\n\n")


