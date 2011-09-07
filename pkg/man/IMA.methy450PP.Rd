\name{IMA.methy450PP}
\alias{IMA.methy450PP}
\alias{methy450batch}
\title{
Data preprocessing and quality control
}
\description{
It allows user to choose several filtering step or modify filtering criteria for specific quality control purpose. These include whether or not to filter probes based on detetion Pvalue; whether or not to remove the loci from the X chromosome; whether or not to transfer the raw \eqn{\beta} vlaue using square arcsine; whether or not to perform normalization; whether or not to remove the loci containing missing \eqn{\beta} values. The user can choose the preprocessing routes and corresponding cutoffs in the argument of this function. 
}
\usage{
IMA.methy450PP(data,na.omit = TRUE,normalization=FALSE,transfm = FALSE,
samplefilterdetectP = c(FALSE,1e-5),samplefilterperc = 0.75,
sitefilterdetectP = c(FALSE,0.05),sitefilterperc = 0.75,
locidiff = c(FALSE,0.01), Xchrom = TRUE)
}

\arguments{
   \item{data}{
an exprmethy450 class returned by the lumi.methy450R function
}
  \item{na.omit}{
if TRUE remove the sites containing missing value
}
  \item{normalization}{
if TRUE, quantile normalization performed
}
  \item{transfm}{
if TRUE, arcsin square root of beta matrix is performed
}
  \item{samplefilterdetectP}{
Default is False, i.e, no sample filtering.Otherwise, chose the cut off of detection Pvalue.
}
  \item{samplefilterperc}{
Remove the sample having over specified percentage of sites with detection
   Pvalue greater than \bold{samplefilterdetectP}.
}
  \item{sitefilterdetectP}{
Default is False, i.e. no sites filtering by detection pvalue. Otherwise,chose the cut off of detection Pvalue.
}
\item{sitefilterperc}{
    Remove the sites having specified percentage of samples with detection
   Pvalue greater than \bold{sitefilterdetectP}.
  }
\item{locidiff}{
if FALSE, don't filter sites by the difference of group \eqn{\beta} value. Otherwise, remove the sites has difference \eqn{\beta} value greater than specited value.
}
  \item{Xchrom}{
if TRUE, Remove the sites on chromosome X
 }
}
\details{
It allows user to choose several filtering step or modify filtering
   criteria for specific quality control purpose. By default, IMA will
   filter out loci with missing \eqn{\beta} value, from the X chromosome
   or with median detection P-value greater than specified value across
   the samples. The option for sample level quality control is also
   provided (B. C. Christensen et al., 2011). Although the raw
   \eqn{\beta} values will be analyzed as recommended by Illumina, the
   user can choose Arcsine squareroot transformation.  Quantile
   normalization is also available for cross sample normalization.
}

\value{
This function will return a methy450batch class including:
\item{bmatrix}{ a matrix of \eqn{\beta} value for individual sites}
\item{detectP}{ a matrix of detection pvalue for individual sites}
\item{annot}{a matrix of annotation information for individual sites}
\item{groupinfo}{ a list of sample ID and penotype of each sample}
\item{TSS1500Ind}{ a list of sites ID belonging to the TSS1500 region of each gene  }
\item{TSS200Ind}{ a list of sites ID belonging to the TSS200 region of each gene }
\item{UTR5Ind}{ a list of sites ID belonging to the 5' UTR region of each gene }
\item{EXON1Ind}{a list of sites ID belonging to the 1st EXON of each gene}
\item{UTR3Ind}{a list of sites ID belonging to the 3' UTR region of each gene}
\item{GENEBODYInd}{a list of sites ID belonging to the gene body region of each gene}
\item{ISLANDInd}{a list of sites ID belonging to the ISLAND region of each UCSC_CPG_ISLAND}
\item{NSHOREInd}{a list of sites ID belonging to the N Shore region of each UCSC_CPG_ISLAND}
\item{SSHOREInd}{a list of sites ID belonging to the S Shore region of each UCSC_CPG_ISLAND}
\item{NSHELFInd}{a list of sites ID belonging to the N Shelf region of each UCSC_CPG_ISLAND}
\item{SSHELFInd}{a list of sites ID belonging to the S Shelf region of each UCSC_CPG_ISLAND}

}
\author{
Dan Wang, Li Yang, Qiang Hu, Dominic J Smiraglia, Song Liu

Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{regionswrapper}},\code{\link{testfunc}},\code{\link{sitetest}},\code{\link{indexregionfunc}}
}
\examples{
setwd(system.file("extdata", package="IMA"))
MethyFileName = "SampleMethFinalReport.txt"
PhenoFileName = "SamplePhenotype.txt"
data = IMA.methy450R(file = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName);
dataf = IMA.methy450PP(data,na.omit = TRUE,normalization=FALSE,transfm =FALSE,samplefilterdetectP = 1e-5,samplefilterperc =0.75,sitefilterdetectP = 0.05,sitefilterperc = 0.75,locidiff = FALSE, Xchrom = TRUE)
#This command will remove the sites have following characters:
#1.containg miss value
#2.with median detection pvalue great than 0.05
#3.on X chromsome
#This command will also remove the samples with 75% sites having detection pvalue greater than 1e-05
}
