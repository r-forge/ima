\name{sitetest}
\alias{sitetest}

\title{
Site level Differential Methylation Analysis
}
\description{
 For each specific site, Wilcoxon rank-sum test (default), Student's t-test ("pooled" or"satterthwaite") and empirical Bayes statistics (G. K. Smyth, 2004) are available for inference in differential testing. Robust linear regression is available as an option to infer methylation change associated with continuous variable (e.g., age). A variety of multiple testing correction algorithms is available, including stringent Bonferroni correction and widely used false discovery rate control. Users can specify the significance criteria in parameter file.  
}
\usage{
sitetest(dataf, gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", 
"limma", "pooled","satterthwaite"), Padj = "BH", concov = "OFF", 
rawpcut = NULL, adjustpcut = NULL, betadiffcut = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
a  methy450batch class returned by the lumi.methy450PP function
}
  \item{gcase}{
 the string name of case group in the sample.txt file
}
  \item{gcontrol}{
the string name of control group in the sample.txt file
}
  \item{testmethod}{
Wilcoxon rank-sum test (default), Student's t-test ("pooled" or "satterthwaite"), or empirical Bayes statistics (G. K. Smyth, 2004).
}
  \item{Padj}{
The methods available for multiple testing correction. Users could choose any method provided in the p.adjust function of R stat package.
}
  \item{concov}{
if "ON", covariates is continuous variable, and the robust linear regression would be used for testing the association between methylation level and continuous phenotye (e.g., age).
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
return a matrix with "P-Value", "Adjust Pval","beta-Difference" in seperate columns for each site
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{regionswrapper}},\code{\link{IMA.methy450PP}},\code{\link{testfunc}},\code{\link{indexregionfunc}}
}
\examples{
#sitetest = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,testmethod = testmethod,Padj=Padj,rawpcut = NULL,adjustpcut =NULL,betadiffcut = NULL)

}
