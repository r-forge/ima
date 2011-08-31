rm(list=ls())
options(stringAsfactors = FALSE)
####Change all the parameter here
####Load data
libPaths = "/home/danwang/R/myR" ####specify your R library
MethyFileName = "SampleMethFinalReport.txt" ###specfiy the original methylation data produced by the GenomeStudio
PhenoFileName = "Phenotype.txt" ###specify the genotype for each sample 
####output file
siteleveltest = "../output/single.test.xls" ###specify the path and name for the site level testing result
#####Preprocessing:lumi.methy450PP()
detectPcut = 0.75 ##Remove the samples with detection Pvalue > 1e-5 across 75% of the sites
filterSample = TRUE ## If this is false, keep all the samples without filtering
na.omit = TRUE ## Remove the sites containing missing beta value
Xchrom = TRUE #### Remove the sites on chromosome X
filterdetectP = TRUE    ## Remove the sites with median detection Pvalue > 0.05 across the samples
normalization = FALSE #### if TRUE, quantile normalization performed
transfm = FALSE ##if TRUE, arcsin square root of beta matrix is performed
locidiff = FALSE ### if TRUE, keep the sites has big variance between two groups.
locidiffcut = 0.01  ### if locidiff = TRUE, keep the sites has variance great than the locidiffcut between two groups 

#####Singtest
testmethod = "limma"  ###other options of testing: limma/pooled/wilcox/satterthwaite for the comparison between two group
concov = "OFF" ### if "ON", covariates is continuous variable
gcase = "g2" ### specify the case group index in the sample.txt file
gcontrol = "g1" ###specify the control group index in the sample.txt file
Padj = "BH" ###specify which method applied for multiple testing correction

###output
outputDES = FALSE
rawpcut = 0.05
adjustpcut = 0.05
betadiffcut = 0.14

####11 region wrapper
sumregion="mean"  ###median/tbrm,specify the way to summarize each region:tbrm(Tukey's Biweight robust average) 
regionint = "TSS1500" ###specify which region is interested in. The options are  eight gene levels("TSS1500","TSS200","UTR5","EXON1","GENEBODY","UTR3") and 5 UCSC CpG Island levels("ISLAND","NSHORE","SSHORE","NSHELF","SSHELF"). 
list11excel = "../output/list11excel.xls"###specify the path and name for region level analysis results
list11Rdata = "../output/list11.Rdata" ###specify the path of Rdata where save the region level analysis results
#####Please install the following R package before you run
#####limma,WriteXLS,preprocessCore,heatmap.2,
######Please download the "methy450package.R" file to your R working directory
.libPaths(libPaths) ##Specify your R library

data = IMA.methy450R(file = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName)#This will take about 5 mins to load the data
#save(data,file="./data/example.raw.data.rda")
#fullannot = list(TSS1500Ind=TSS1500Ind,TSS200Ind=TSS200Ind,UTR5Ind=UTR5Ind,EXON1=EXON1Ind, GENEBODYInd=GENEBODYInd, UTR3Ind=UTR3Ind,ISLANDInd=ISLANDInd,NSHOREInd=NSHOREInd,SSHOREInd=SSHOREInd,NSHELFInd=NSHELFInd,SSHELFInd=SSHELFInd)
#save(fullannot,file="./data/fullannotInd.rda")
#fullannot[["TSS200Ind"]][1000:1010]
dataf = IMA.methy450PP(data,na.omit = na.omit,normalization=normalization,transfm = transfm,filterdetectP = filterdetectP, detectPcut = detectPcut, Xchrom = Xchrom)  ##About 2 mins
slotNames(dataf)
dataf@TSS200Ind[1000:1010]
#save(dataf,file="./data/example.filtered.data.rda")
############Single probe test with the "BH" adjustment
singtest1 = sigtest(dataf,gcase=gcase,gcontrol=gcontrol,test = testmethod,Padj=Padj,outputDES = FALSE,rawpcut = NULL,adjustpcut =NULL,betadiffcut = NULL)
singtest = outputDESfunc(singtest1,outputDES = TRUE,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14)
DESsingtest = as.data.frame(singtest1);
write.table(DESsingtest,file=siteleveltest,row.names=TRUE)
###WriteXLS not working on the data exceeds 65535 rows or 256 columns
###########11 region wrapper
regionswrapper(dataf,sumregion=sumregion,gcase=gcase,gcontrol=gcontrol,testmethod=testmethod,Padj=Padj,concov=concov,list11excel=list11excel,list11Rdata =list11Rdata,outputDES = TRUE,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14)###This one will take one hour
###Interested in one specific region level analysis: e.g. TSS1500
betar = sumregionfun(indexlist=dataf@TSS1500Ind[1:100],beta=dataf@bmatrix,sumregion=sumregion)
group = dataf@groupinfo
grouplev = group$group[match(colnames(betar),group$samplename)]
print(grouplev)
TSS1500test = testfunc(eset = betar,testmethod=testmethod,Padj=Padj,concov=concov,grouplev = grouplev,gcase = "g2",gcontrol="g1")
TSS1500DES = outputDESfunc(TSS1500test,outputDES = TRUE,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14)

####Zoom in the methylation on individual genes/cpG Islands 
####i.e. 
load("./data/fullannotInd.rda")##this file could download from the R-forge.
indlists = c("BRCA1", "MLH1", "CCNE1", "PTEN", "PALB2")
index = match("TSS1500Ind",names(fullannot))
annot = fullannot[[index]]
index = match(indlists,names(annot))
indexlist = annot[index]     
beta = data@bmatrix;
eset = sumregionfun(indexlist,beta,"mean");
testfunc(eset,concov = "OFF",testmethod="limma",Padj="BH",gcase ="g2",gcontrol="g1",grouplev=grouplev)


