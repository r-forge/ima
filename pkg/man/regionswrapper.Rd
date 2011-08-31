\name{regionswrapper}
\alias{regionswrapper}
\title{
Differentially test on the 11 annotated region-level analysis separately
}
\description{
this function will output all the statistical results for the 11
annotated region-level analysis result in separated sheet of excel
file, and the result could be used for further analysis such as GO
functional annotation, pathway analysis, network etc. This function will
took up to 2 hours to run.
}
\usage{
regionswrapper(dataf, sumregion = c("mean", "median", "tbrm"), gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = "BH", concov = c("OFF", "ON"), list11excel = list11excel, list11Rdata = list11Rdata, outputDES = FALSE, rawpcut = NULL, adjustpcut = NULL, betadiffcut = NULL)
}

\arguments{
  \item{dataf}{
a  methy450batch class returned by the lumi.methy450PP function
}
  \item{sumregion}{
the methods can be used to summarize the targed sites for each region
}
  \item{gcase}{
the string name of case group in the sample.txt file
}
  \item{gcontrol}{
the string name of control group in the sample.txt file
}
  \item{testmethod}{
Wilcoxon rank-sum test (default), Student's t-test and empirical Bayes statistics (G. K. Smyth, 2004).
}
  \item{Padj}{
the method applied for multiple testing correction.users could choose any method in the p.adjust function such as FDR, BH etc.
}
  \item{concov}{
 if "ON", covariates is continuous variable, then the robust linear regression would be used for testing. 
}
  \item{list11excel}{
the name of output excel file containing the differential testing result for each annotated region
}
  \item{list11Rdata}{
the name of Rdata containing the differential testing result for each annotated region
}
 \item{outputDES}{
if outputDES = TRUE, only out the differential methylated sites by using
defined cut off
}
  \item{rawpcut}{
only used when outputDES = TRUE, and define the cut off of the raw pvalue
}
\item{adjustpcut}{
 only used when outputDES = TRUE, and define the cut off of the adjusted pvalue

}
  \item{betadiffcut}{
only used when outputDES = TRUE, and define the cut off of the beta
  difference between two groups.
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
  \code{\link{IMA.methy450PP}},\code{\link{IMA.methy450R}},\code{\link{testfunc}},\code{\link{sigtest}},\code{\link{sumregionfun}}
}
\examples{
#regionswrapper(dataf,sumregion=sumregion,gcase=gcase,gcontrol=gcontrol,
#testmethod=testmethod,Padj=Padj,concov=concov,list14excel="../output/list14result.xls")
}