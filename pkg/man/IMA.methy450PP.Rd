\name{IMA.methy450PP}
\alias{IMA.methy450PP}
\alias{methy450batch}
\title{
Data preprocessing and quality control
}
\description{
It allows user to choose several filtering step or modify filtering criteria for specific quality control purpose. These include whether or not to filter probes based on detetion Pvalue; whether or not to remove the loci from the X chromosome; whether or not to transfer the raw \beta vlaue using square arcsine; whether or not to perform normalization; whether or not to remove the loci containing missing \beta values. The user can choose the preprocessing routes and corresponding cutoffs in the argument of this function. 
}
\usage{
IMA.methy450PP(data, na.omit = TRUE, normalization = FALSE, transfm = FALSE, filtersample = TRUE, filterdetectP = TRUE, detectPcut = 0.75, locidiff = FALSE, locidiffcut = 0.01, Xchrom = TRUE)
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
  \item{filtersample}{
Keep the samples with detection Pvalue <10e-5 on more than 75\% sites 
}
  \item{filterdetectP}{
Remove the sites with detection Pvalue > 0.01 across  detectPcut percentage of the samples
}
  \item{detectPcut}{
Remove the sites with detection Pvalue > 0.01 across 75\% of the samples
}
\item{locidiff}{
if TRUE, keep the sites has big variance between two groups.
}
  \item{locidiffcut}{
if locidiff = TRUE, keep the sites has variance great than the locidiffcut between two groups
}
  \item{Xchrom}{
if TRUE, Remove the sites on chromosome X
 }
}
\details{
It allows user to choose several filtering step or modify filtering criteria for specific quality control purpose. By default, IMA will filter out loci with missing beta value, from the X chromosome or with median detection P-value greater than 0.05 across the samples. The option for sample level quality control is also provided (B. C. Christensen et al., 2011). Although the raw beta values will be analyzed as recommended by Illumina, the user can choose Arcsine squareroot transformation.  Both scaling and quantile normalization options are available for cross sample normalization.
}

\value{
This function will return a methy450batch class including:
\item{bmatrix}{ a matrix of beta value for individual sites}
\item{detectP}{ a matrix of detection pvalue for individual sites}
\item{annot}{a matrix of annotation information for individual sites}
\item{groupinfo}{ a list of sample ID and penotype of each sample}
\item{TSS1500Ind}{ a list of sites ID belonging to the TSS1500 region of each gene  }
\item{TSS200Ind}{ a list of sites ID belonging to the TSS200 region of each gene }
\item{UTR5Ind}{ a list of sites ID belonging to the 5' UTR region of each gene }
\item{EXON1Ind}{a list of sites ID belonging to the 1st EXON of each gene}
\item{UTR3Ind}{a list of sites ID belonging to the 3' UTR region of each gene}
\item{GENEBODYInd}{a list of sites ID belonging to the gene body region of each gene}
\item{ISLANDInd}{a list of istes ID belonging to the ISLAND region of each UCSC_CPG_ISLAND}
\item{ NSHOREInd}{a list of istes ID belonging to the N Shore region of each UCSC_CPG_ISLAND}
\item{SSHOREInd}{a list of istes ID belonging to the S Shore region of each UCSC_CPG_ISLAND}
\item{NSHELFInd}{a list of istes ID belonging to the N Shelf region of each UCSC_CPG_ISLAND}
\item{SSHELFInd}{a list of istes ID belonging to the S Shelf region of each UCSC_CPG_ISLAND}

}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{regionswrapper}},\code{\link{testfunc}},\code{\link{sigtest}},\code{\link{sumregionfun}}
}
\examples{
#dataf =IMA.methy450PP(data,na.omit=na.omit,normalization=normalization,
#transfm = transfm,filterdetectP = filterdetectP, detectPcut =detectPcut, Xchrom = xchrom)
#About 2 mins for preprocessing ~450k sites
}
