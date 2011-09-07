\name{outputDMfunc}
\alias{outputDMfunc}
\title{
Output the differentially methylated sites/regions
}
\description{
  This function will output the testing result with user-specified significance criteria. 
 }
\usage{
outputDMfunc(out, rawpcut = 0.05, adjustpcut = 0.05, betadiffcut = 0.14)
}
\arguments{
  \item{out}{
this is the output by function of sitetest, testfunc or regionwrapper with three columns for each site/region: pvalue, adjusted pvalue and
\eqn{\beta} difference.
}
  \item{rawpcut}{
either null or define the cut off of the raw pvalue
}
\item{adjustpcut}{
either null or define the cut off of the adjusted pvalue

}
  \item{betadiffcut}{
either null or define the cut off of the \eqn{\beta} difference between two groups.
}

}

\value{
return a matrix with "P-Value", "Adjust Pval","beta-Difference" in seperate columns for each site/region
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{regionswrapper}},\code{\link{IMA.methy450PP}},\code{\link{testfunc}},\code{\link{indexregionfunc}}
}
\examples{
#sitetest = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,test = testmethod,Padj=Padj,outputDES = FALSE,rawpcut = NULL,adjustpcut =NULL,betadiffcut = NULL)
#sitetest = outputDMfunc(sitetest,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14)
}
