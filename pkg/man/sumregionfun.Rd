\name{sumregionfun}
\alias{sumregionfun}
\title{
  Methylation Index Calculation
}
\description{
  for each specific region (e.g., promoter), IMA will collect the loci within it and derive an index of overall region methylation value. Currently, there are three different index metrics implemented in IMA: mean, median, and Tukey's Biweight robust average.  By default, the median beta values will be used as the region's methylation index for further analysis. 
}
\usage{
sumregionfun(indexlist = dataf@GENEInd, beta, sumregion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexlist}{specify which annotated region you are interested
}
\item{beta}{
  beta matrix for each ~450k site
}
\item{sumregion}{
  The method derive the overall region methylation value for each gene or UCSC_CPG_ISLAND. The options are mean, median, and Tukey's Biweight robust average.
}
}
\value{
  return a summarized matrix for the annotated region
}
\references{
  https://r-forge.r-project.org/projects/ima/
}
\author{
  Dan Wang; Li Yan; Qiang Hu;Dominic J Smiraglia; Song Liu
  Maintainer: Dan Wang <wangdan412@gmail.com>
}
\examples{
#data(fullannot)
#indlists = c("BRCA1", "MLH1", "CCNE1", "PTEN", "PALB2")
#index = match("TSS1500Ind",names(fullannot))
#annot = fullannot[[index]]
#index = match(indlists,names(annot))
#indexlist = annot[index]     
#beta = data@bmatrix;
#eset = sumregionfun(indexlist,beta,"mean");
#testfunc(eset,concov = "OFF",testmethod="limma",Padj="BH",gcase ="g2",gcontrol="g1",grouplev=grouplev)
}