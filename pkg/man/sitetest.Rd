\name{sitetest}
\alias{sitetest}

\title{
Differential Methylation Analysis
}
\description{
 For each specific region, Wilcoxon rank-sum test (default), Student's t-test and empirical Bayes statistics (G. K. Smyth, 2004) are available for inference in differential testing. Robust linear regression is available as an option to infer methylation change associated with continuous variable (e.g., age). A variety of multiple testing correction algorithms is available, including stringent Bonferroni correction and widely used false discovery rate control. Users can specify the significance criteria in parameter file. The same statistical inference and multiple test correction procedures described above will also be applied to each single site to obtain site-level differential methylation inference.  
}
\usage{
sitetest(dataf, gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = "BH", concov = "OFF", rawpcut = NULL, adjustpcut = NULL, betadiffcut = NULL)
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
Wilcoxon rank-sum test (default), Student's t-test and empirical Bayes statistics (G. K. Smyth, 2004).
}
  \item{Padj}{
 specify which method applied for multiple testing correction.
}
  \item{concov}{
if "ON", covariates is continuous variable
}
  \item{rawpcut}{
define the cut off of the raw pvalue
}
\item{adjustpcut}{
define the cut off of the adjusted pvalue

}
  \item{betadiffcut}{
define the cut off of the beta difference between two groups.
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
#sitetest = outputDESfunc(sitetest,outputDES = TRUE,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14)
}
