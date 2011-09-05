\name{outputDESfunc}
\alias{outputDESfunc}
\title{
Output the differential expressed sites/regions
}
\description{
  This function will subset the testing result by using the user
  specified cut off. 
 }
\usage{
outputDESfunc(out, rawpcut = 0.05, adjustpcut = 0.05, betadiffcut = 0.14)
}
\arguments{
  \item{out}{
this is the return results by function of sigtest,
testfunc,regionwrapper with three columns: pvalue,adjusted pvalue and
beta difference.
}
  \item{rawpcut}{
define the cut off of the raw pvalue
}
\item{adjustpcut}{
define the cut off of the adjusted pvalue

}
  \item{betadiffcut}{
define the cut off of the \eqn{\beta} difference between two groups.
}

}

\value{
return a matrix with "P-Value", "Adjust Pval","beta-Difference" in columns for each site
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{regionswrapper}},\code{\link{IMA.methy450PP}},\code{\link{testfunc}},\code{\link{sumregionfun}}
}
\examples{
#sitetest = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,test = testmethod,Padj=Padj,outputDES = FALSE,rawpcut = NULL,adjustpcut =NULL,betadiffcut = NULL)
#sitetest = outputDESfunc(sitetest,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14)
}
