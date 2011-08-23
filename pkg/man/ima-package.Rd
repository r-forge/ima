\name{IMA}
\alias{Illumina Methylation Analyzer}
\docType{package}
\title{
IMA(Illumina Methylation Analyzer) is a package designed to auto-mate the pipeline for analyzing site-level and region-level methylation changes in epigenetic studies utilizing the 450K DNA methylation microarray
}
\description{
IMA is implemented in R and can be run on any platform with an existing R and Bioconductor installation. The user can either run the pipeline with default setting or specify optional routes in the parameter file. IMA takes as input the beta values representing the methylation levels of individual site reported by Illumina BeadStudio or GenomeStudio software. It allows user to choose several filtering step or modify filtering criteria for specific quality control purpose.For each specific region (e.g., promoter), IMA will collect the loci within it and derive an index of overall region methylation value.Then those value could be used to do differential testing by different method.
}
\details{
	\tabular{ll}{
	Package: \tab IMA\cr
	Type: \tab Package\cr
	Version: \tab 1.0\cr
	Date: \tab 2011-08-22\cr
	License: \tab GPL\cr
	LazyLoad: \tab yes\cr
	}
}
\author{
Dan Wang, Li Yang, Qiang Hu,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\references{
  https://r-forge.r-project.org/projects/ima/
}
\keyword{ package }
\examples{
#data = lumi.methy450R(file = inputmethydata,groupfile = genotypeinput)
#dataf = lumi.methy450PP(data,na.omit =na.omit,normalization=normalization,transfm = transfm,filterdetectP =filterdetectP, detectPcut = detectPcut, Xchrom = xchrom)
#singtest1 = sigtest(dataf,gcase=gcase,gcontrol=gcontrol,test =testmethod,Padj=Padj)
#regionswrapper(dataf,sumregion=sumregion,gcase=gcase,gcontrol=gcontrol,testmethod=testmethod,Padj=Padj,concov=concov,list14excel="../output/list14result.xls")
}
