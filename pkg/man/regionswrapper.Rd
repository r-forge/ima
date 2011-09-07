\name{regionswrapper}
\alias{regionswrapper}
\title{
Differential methylation test on all 11 annotated regions
}
\description{
This function will output all the differential methylation test results for the 11
annotated regions, and store the result in separated sheets of an excel
file.
}
\usage{
regionswrapper(dataf, indexmethod = c("mean", "median", "tbrm"), gcase = "g2", 
gcontrol = "g1",testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), 
Padj = "BH", concov = c("OFF", "ON"),list11excel,list11Rdata, rawpcut = NULL,
adjustpcut = NULL, betadiffcut = NULL)
}

\arguments{
  \item{dataf}{
a  methy450batch class returned by the lumi.methy450PP function
}
  \item{indexmethod}{
the methods available to derive an index of overall methylation level for each region
}
  \item{gcase}{
the string name of case group in the sample.txt file
}
  \item{gcontrol}{
the string name of control group in the sample.txt file
}
  \item{testmethod}{
Wilcoxon rank-sum test (default), Student's t-test ("pooled" or "satterthwaite") or empirical Bayes statistics (G. K. Smyth, 2004).
}
  \item{Padj}{
the method available for multiple testing correction. Users could choose any method provided in the p.adjust function of R stat package.
}
  \item{concov}{
 if "ON", covariates is continuous variable, and the robust linear regression would be used for testing the association between methylation level and continuous phenotye (e.g., age). 
}
  \item{list11excel}{
the name of output excel file containing the differential testing result for each annotated region
}
  \item{list11Rdata}{
the name of Rdata containing the differential testing result for each annotated region
}
  \item{rawpcut}{
either null or define the cut off of the raw pvalue
}
\item{adjustpcut}{
either null or define the cut off of the adjusted pvalue

}
  \item{betadiffcut}{
either null or define the cut off of the beta difference between two groups.
}
 
}


\value{
Return an excel file containing the differential testing result for each annotated region and a Rdata containing the differential testing  result for each annotated region
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu

Maintainer: Dan Wang <wangdan412@gmail.com>
}

\seealso{
  \code{\link{IMA.methy450PP}},\code{\link{IMA.methy450R}},\code{\link{testfunc}},\code{\link{sitetest}},\code{\link{indexregionfunc}}
}
\examples{
setwd(system.file("extdata", package="IMA"))
MethyFileName = "SampleMethFinalReport.txt"
PhenoFileName = "SamplePhenotype.txt"
data = IMA.methy450R(file = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName);
dataf = IMA.methy450PP(data,na.omit = TRUE,normalization=FALSE,transfm = FALSE,samplefilterdetectP = 1e-5,samplefilterperc = 0.75,sitefilterdetectP = 0.05,sitefilterperc = 0.75,locidiff = FALSE, Xchrom = TRUE)
regionswrapper(dataf,indexmethod ="mean",gcase = "g2",gcontrol="g1",testmethod = "wilcox",Padj="BH",concov = "OFF",list11excel="list14result.xls",list11Rdata="list14result.Rdata",rawpcut = NULL,adjustpcut = NULL,betadiffcut = NULL)
}
