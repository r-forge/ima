\name{IMA.methy450R}
\alias{IMA.methy450R}
\alias{exprmethy450}
\title{
Read in methylation 450k data
}
\description{
 This function read in methylation 450k data by one single command. Twoinput data required. One is the methylation expression .txt data which
 exported from BeadStudio with tab separated.The other one is the phenotype data which prepared by the users. This function also return a cluster tree to show the unsuperwised clustering of samples by
using all 450k sites.The Quality control figures for each samples are also returned in the QC.pdf file and can be used to decide sample/sites filtering.
}
\usage{
IMA.methy450R(file = MethtFileName, columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"), groupfile = PhenoFileName)
}

\arguments{
  \item{file}{
This is the file paths/name containing 450k methylation sites with
correspinding beta value,detection pvalue and annotation information and was produced from the illumina BeadStudio or GenomeStudio
software.
}
  \item{columnGrepPattern}{
specify the string pattern's name of beta value and detection pvalue in the .txt file
}
  \item{groupfile}{
This is the file paths/name of the penotype data.
}
}

\value{
This function will return a exprmethy450 class including below values:a beta matrix, a annot matrix, a matrix with detect_pvalue for each ~450k site and a group info
matrix.
\item{bmatrix}{ beta matrix for each ~450k site on each sample}
\item{detectP}{ a matrix of detection pvalue for each ~450k site on each sample}
\item{annot}{a matrix of annotation for each ~450k site}
\item{groupinfo}{ a list of sample ID and penotype for each sample}
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450PP}},\code{\link{regionswrapper}},\code{\link{testfunc}},\code{\link{sigtest}},\code{\link{sumregionfun}}
}
\examples{
setwd(system.file("extdata", package="IMA"))
MethyFileName = "SampleMethFinalReport.txt"
PhenoFileName = "Phenotype.txt"
data = IMA.methy450R(file = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName);
}
