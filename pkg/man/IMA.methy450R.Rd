\name{IMA.methy450R}
\alias{IMA.methy450R}
\alias{exprmethy450}
\title{
Load methylation 450k data
}
\description{
This function will load the methylation 450k data. The input for the package consists basically of two files containing  value methylation data (including annotation) produced by BeadStudio or GenomeStudio, and sample phenotype data prepared by the user. \emph{IMA.methy450R} load the input files with a single command described below and a \emph{exprmethy450} object will be created, which includes the following features: \eqn{\beta} value matrix, locus annotion, detection Pvalue, sample phenotype information. Besides, baisc quality control information will be outputed in the QC.pdf, which include unsupervised sample clustering using all loci, boxplot for  value of each sample, and barplot showing the percent of loci with detection Pvalue smaller than 1e-5 in each sample.
}
\usage{
IMA.methy450R(file , columnGrepPattern = list(beta =".AVG_Beta"
, detectp = ".Detection.Pval"), groupfile)
}

\arguments{
  \item{file}{
This is the input file containing 450k methylation data including \eqn{\beta} value, detection pvalue and annotation information. The file can be be produced from the illumina BeadStudio or GenomeStudio
software.
}
  \item{columnGrepPattern}{
Specify the columns with corresponding characters from the file produced by llumina BeadStudio or GenomeStudio software will be used. 
}
  \item{groupfile}{
This is the input file containing the penotype data.
}
}

\value{
This function will return an \emph{exprmethy450} class including these
features: a \eqn{\beta} value matrix, an annotation matrix, a detection
pvalue matrix and a list of sample phenotype. A QC.pdf file will be
output as well which including the basic quality control infromation.
\item{bmatrix}{ \eqn{\beta} value matrix for each site in each sample}
\item{detectP}{ a matrix of detection pvalue for each ~450k site on each sample}
\item{annot}{a matrix of annotation for each ~450k site}
\item{groupinfo}{ a list of sample ID and penotype for each sample}
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
  
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450PP}},\code{\link{regionswrapper}},\code{\link{testfunc}},\code{\link{sitetest}},\code{\link{indexregionfunc}}
}
\examples{
setwd(system.file("extdata", package="IMA"))
MethyFileName = "SampleMethFinalReport.txt"
PhenoFileName = "SamplePhenotype.txt"
data = IMA.methy450R(file = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName);
}
