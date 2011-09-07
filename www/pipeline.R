rm(list=ls())
options(stringAsfactors = FALSE)

#############################################################################################################
#############################################################################################################
########################Prepare all the parameters here#######################################################
#############################################################################################################
#############################################################################################################

#######################options in IMA.methy450R#############################################################
######################Load data#############################################################################
libPaths = "/home/danwang/R/myR" ####specify the location of your R library
MethyFileName = "SampleMethFinalReport.txt" ###specfiy the original methylation data produced by the GenomeStudio
PhenoFileName = "SamplePhenotype.txt" ###specify the phenotype for each sample
#############################################################################################################

###########output file#####################################################################################
siteleveltest = "./sitelevle.test.xls" ###specify the path and name for the site level testing result
list11excel = "../list11excel.xls"###specify the path and name for region level analysis results
list11Rdata = "../list11.Rdata"###specify the path of Rdata which stores the region level analysis results
############################################################################################################

#################Preprocessing:IMA.methy450PP ##############################################################
samplefilterdetectP = 1e-5 ## The cutoff for sample level detection Pvalue
samplefilterperc = 0.75    ## The percent of loci with detection Pvalue less than samplefilterdetectP in each sample
sitefilterdetectP = 0.05   ## The cutoff for site level detection Pvalue
sitefilterperc = 0.5       ## The percent of sample with detection Pvalue less than sitefilterdetectP for each site
na.omit = TRUE ## Remove the sites containing missing beta value
Xchrom = TRUE #### Remove the sites on chromosome X
normalization = FALSE #### if TRUE, quantile normalization performed
transfm = FALSE ##if TRUE, arcsin square root of beta matrix is performed
locidiff = FALSE ### if not FALSE, keep the sites with specified beta value difference.
################################################################################################################

############sitetest/regionwrapper####################################################################################
testmethod = "limma"  ###other options of testing: wilcox/pooled/satterthwaite for the comparison between two group
concov = "OFF" ### if "ON", covariates is continuous variable
gcase = "g2" ### specify the case group index in the sample.txt file
gcontrol = "g1" ###specify the control group index in the sample.txt file
Padj = "BH" ###specify which method applied for multiple testing correction. The user can choose the methods provided by p.adjust function of R stat package 
indexmethod ="mean"  ###median/tbrm, specify the way to derive an index of overall methylation value of each region. "tbrm" is Tukey's Biweight robust average. 
###################################################################################################################

####################################output the differential test######################################################
rawpcut = 0.05     ## cut off for raw pvalue
adjustpcut = 0.05  ## cut off for adjusted p value
betadiffcut = 0.14 ## cut off for beta value difference
######################################################################################################################

######################################################################################################################
###############################End of the parameter specification###########################################################
######################################################################################################################


################################ Analysis Routes ####################################################################
######################################################################################################################

.libPaths(libPaths) ##Specify your R library
data = IMA.methy450R(file = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName)# load the data
dataf = IMA.methy450PP(data,na.omit = TRUE,normalization=FALSE,transfm = FALSE,samplefilterdetectP = 1e-5,samplefilterperc = 0.75,sitefilterdetectP = 0.05,locidiff = FALSE, Xchrom = TRUE)## QC filtering

sitetest = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,testmethod = testmethod,Padj=Padj,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut) ###site level test with the "BH" adjustment
write.table(sitetest,file=siteleveltest,row.names=TRUE) # saving the reults (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)

regionswrapper(dataf,indexmethod =indexmethod,gcase = gcase,gcontrol=gcontrol,testmethod = testmethod,Padj=Padj,concov = concov,list11excel=list11excel,list11Rdata=list11Rdata,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut) # region level testing for all 11 categories of annotated regions
######################################################################################################################
######################################################################################################################


