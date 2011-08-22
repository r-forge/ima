rm(list=ls())
####Change all the parameter here
####Load data
libPaths = "/home/danwang/R/myR" ####specify your R library
sourcefile = "./methy450package.R"  ###specify the path and file where saved the "methy450package.R"
inputmethydata = "../data/SampleMethFinalReport_Smiraglia.txt" ###specfiy the original methydata produced by the GenomeStudio
genotypeinput = "../data/sample.txt" ###specify the genotype for each sample 
####output file
siteleveltest = "../output/single.test.xls" ###specify the path and name for the site level testing
#####Preprocessing:lumi.methy450PP()
filterSample = TRUE ## keep the samples with detection Pvalue <10e-5 on more than 75% sites  
na.omit = TRUE ## Remove the sites containing missing beta value
xchrom = TRUE #### Remove the sites on chromosome X
filterdetectP = TRUE    ## Remove the sites with detection Pvalue > 0.01 across  detectPcut percentage of the samples
detectPcut = 0.75 ##Remove the sites with detection Pvalue > 0.01 across 75% of the samples
normalization = FALSE #### if TRUE, quantile normalization performed
transfm = FALSE ##if TRUE, arcsin square root of beta matrix is performed
locidiff = FALSE ### if TRUE, keep the sites has big variance between two groups.
locidiffcut = 0.01  ### if locidiff = TRUE, keep the sites has variance great than the locidiffcut between two groups 

#####Singtest
testmethod = "wilcox"  ###other options of testing: limma/pooled/satterthwaite for the comparison between two group
concov = "OFF" ### if "ON", covariates is continuous variable
gcase = "g2" ### specify the case group index in the sample.txt file
gcontrol = "g1" ###specify the control group index in the sample.txt file
Padj = "BH" ###specify which method applied for multiple testing correction
####14 region wrapper
sumregion="mean"  ###median/tbrm,specify the way to summarize each region:tbrm(Tukey's Biweight robust average) 
regionint = "PROMOTER" ###specify which region is interested in

#####Please install the following R package before you run
#####limma,WriteXLS,preprocessCore,heatmap.2,
######Please download the "methy450package.R" file to your R working directory
.libPaths(libPaths) ##Specify your R library
source(sourcefile)
data = lumi.methy450R(file = inputmethydata,groupfile = genotypeinput)#This will take about 10 mins to load the data 
save(data,file = "../Rdata/raw.Rdata")
dataf = lumi.methy450PP(data,na.omit = na.omite,normalization=normalization,transfm = transfm,filterdetectP = filterdetectP, detectPcut = detectPcut, Xchrom = xchrom)  ##About 2 mins
############Single probe test with the "BH" adjustment
singtest1 = sigtest(dataf,gcase=gcase,gcontrol=gcontrol,test = testmethod,Padj=Padj)
singtest1 = as.data.frame(singtest1);
write.table(singtest1,file=siteleveltest,row.names=TRUE)
###WriteXLS not working on the data exceeds 65535 rows or 256 columns
###########14 region wrapper
regionswrapper(dataf,sumregion=sumregion,gcase=gcase,gcontrol=gcontrol,testmethod=testmethod,Padj=Padj,concov=concov,list14excel="../output/list14result.xls")###This one will take one hour

###Interested in annoted region level: e.g. PROMOTER
eval(parse(text = paste("indexlist = dataf@",regionint,"Ind",sep="")))
betar = sumregionfun(indexlist=indexlist,beta=dataf@bmatrix,sumregion=sumregion)
group = dataf@groupinfo
grouplev = group$group[match(colnames(beta),group$samplename)]
PROMOTERtest = testfunc(eset = betar,sumregion=sumregion,testmethod=testmethod,Padj=Padj,concov=concov,lev1 = which(grouplev==gcase),lev2 = which(grouplev==gcontrol),grouplev = grouplev)

