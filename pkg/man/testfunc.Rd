\name{testfunc}
\alias{testfunc}
\title{
  Differential methylation testing function
}
\description{
  Choose the testing function for differential methylation inference 
}

\usage{
testfunc(eset, concov = c("ON", "OFF"), testmethod = c("wilcox", "limma", 
"pooled", "satterthwaite"),Padj = c("holm", "hochberg", "hommel", "bonferroni", 
"BH", "BY", "fdr", "none"), grouplev, gcase = "g1",gcontrol = "g2")
}

\arguments{
  \item{eset}{
a beta matrix
}
\item{concov}{
  if "ON", covariates is continuous variable,"grouplev" should be the
  continuous observation. 

}
  \item{testmethod}{
Wilcoxon rank-sum test (default), Student's t-test ("pooled" or "satterthwaite"), or empirical Bayes statistics (G. K. Smyth, 2004).
}
  \item{Padj}{
The methods available for multiple testing correction. Users could choose any method provided in the p.adjust function of R stat package.
}
  \item{grouplev}{
specify the vector for covariate:if concov is "ON", this should be the
  continuous values, if concov is "OFF", this should be categorical values. 
}
\item{gcase}{
the string name of case group in the sample.txt file
}
\item{gcontrol}{
the string name of control group in the sample.txt file
}

}
\value{
return a matrix with "P-Value", "Adjust Pval","beta-Difference" in seperate columns for each site or annotated region
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu

Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{regionswrapper}},\code{\link{IMA.methy450PP}},\code{\link{sitetest}},\code{\link{indexregionfunc}}
}
\examples{
#data("fullannotInd")
#indlists = c("BRCA1", "MLH1", "CCNE1", "PTEN", "PALB2")
#annot = fullannot[[match("TSS1500Ind",names(fullannot))]]
#indexlist = annot[match(indlists,names(annot))]     
#eset = indexregionfunc(indexlist,data@bmatrix,"mean");
#testfunc(eset,concov = "OFF",testmethod="limma",Padj="BH",gcase ="g2",gcontrol="g1",
#grouplev=grouplev)
}
