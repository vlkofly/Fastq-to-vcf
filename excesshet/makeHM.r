# create excess heterozygosity mask
# you need to change individual names and population names
setwd("/home/aa/alpine/halleriGenome/makeHM")
# make a name variable modified by Jakub 2021
library(data.table,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
names<-c("CHROM","POS","GRO_01dl","GRO_02dl","GRO_03dl","GRO_04dl","GRO_05dl","KAR_01dl","KAR_02dl","KAR_03dl","KAR_04dl","KAR_05dl","KAR_06dl","LOM_01dl","LOM_02dl","LOM_03dl","LOM_04dl","LOM_05dl","LOM_06dl","LOM_07dl","LOM_08dl","LOM_09dl","MAL_01dl","MAL_02dl","MAL_03dl","MAL_04dl","MAY_01dl","MAY_02dl","MAY_03dl","MAY_04dl","MAY_05dl","MAY_06dl","NKM_01dl","NKM_02dl","NKM_03dl","NKM_04dl","NKM_05dl","NKM_06dl","NKM_07dl","NKM_08dl","NKM_09dl","OSL_02dl","OSL_03dl","OSL_04dl","OSL_05dl","PIZ_02dl","PIZ_03dl","PIZ_04dl","PIZ_05dl","PIZ_06dl","PIZ_08dl","PIZ_09dl","PIZ_11dl","PLE_01dl","PLE_02dl","PLE_03dl","PLE_04dl","PLE_05dl","PLE_06dl","PLE_07dl","PLE_08dl","PLE_09dl","PLE_10dl","PLE_11dl","POT_01dl","POT_02dl","POT_03dl","POT_04dl","SPN_01dl","SPN_02dl","SPN_03dl","SPN_04dl","SPN_05dl","SPN_06dl","SPN_07dl","SPN_08dl","SPN_09dl","SPN_10dl","SPN_11dl","SPN_12dl","SPN_13dl","STD_01dl","STD_02dl","STD_03dl","STD_04dl","STD_05dl","STD_06dl","STD_07dl","STD_08dl","STR_01dl","STR_02dl","STR_03dl","STR_04dl","STR_05dl","STR_06dl","STR_07dl","STR_08dl","STR_09dl","STU_01dl","STU_02dl","STU_03dl","STU_04dl","STU_05dl","STU_06dl","SUN_01dl","SUN_02dl","SUN_03dl","SUN_04dl","SUN_05dl","VLH_01dl","VLH_02dl","VLH_03dl","VLH_05dl","VLH_06dl","VLH_07dl","VLH_08dl","VLH_09dl","VLH_10dl","VLH_11dl","VLH_12dl")
a<-fread("bibp.table",col.names = names )
pop<-c("GRO","KAR","LOM","MAL","MAY","NKM","OSL","PIZ","PLE","POT","SPN","STD","STR","STU","SUN","VLH")
popsum<-function(p){
	colind<-grep(p,colnames(a)) # get column index of particular population
	numind<-length(colind) # get number of individuals per population
	nameind<-colnames(a)[grep(p,colnames(a))] # name of columns
	pt<-a[,..nameind] # subset the data.table
	pt<-cbind(a[,c(1,2)],pt) # add CHROM and POS and GENE
	pt$het<-apply(pt=="0/1",1,sum) # how many hets per site
	pt$mis<-apply(pt=="./.",1,sum) # how many missing per site
	pt$hm<-(numind-pt$mis-pt$het) # number of homozygots
	pt$pmis<-pt$mis/numind # percentage of missingness
	pt<-subset(pt,pt$pmis<0.25) # remove sites with high missingness
	pout<-pt[,c("CHROM","POS","hm","pmis")] # then we will work only with the number of homozygots and percentage of missingness
	colnames(pout)[c(3,4)]<-c(paste(p,"hm",sep=""),paste(p,"pmis",sep=""))
	return(pout)

}

# now run this function per population
popres<-lapply(pop,popsum)
#merge the result into one table
poptable<-Reduce(merge,popres)
head(poptable)
dim(poptable)
sapply(popres,dim)
idhom<-grep("hm",colnames(poptable)) # get indices of homozygosity columns
hm<-poptable[,..idhom]
poptable$fixhet<-apply(hm==0,1,sum)
fixhetpersite<-subset(poptable,poptable$fixhet>0)
fixhetpersite<-fixhetpersite[,c("CHROM","POS","fixhet")]
write.table(fixhetpersite,"fixedHetPerSite.txt",quote=F,row.names=F)
r<-read.table("genesLyV2.gff",h=F)# genes
r$fixhet<-apply(r,1,function(x) nrow(subset(fixhetpersite,fixhetpersite$CHROM==x[1] & fixhetpersite$POS %in% x[2]:x[3])))
write.table(r,"fixedHetPerGene.txt",quote = F,row.names = F)
a5<-subset(fixhetpersite,fixhetpersite$fixhet > 1)
r$fixhet<-apply(r,1,function(x) nrow(subset(a5,a5$CHROM==x[1] & a5$POS %in% x[2]:x[3])))
write.table(r,"fixedHetPerGene_2orMorePop.txt",quote = F,row.names = F)
#and 5 or more positions 
r1<-subset(r,r$fixhet>4)
write.table(r1,"fixedHetPerGene_2orMorePop_5orMoreSites.txt",quote = F,row.names = F,col.names = F)



