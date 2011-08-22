\name{lumi.methy450PP}
\alias{lumi.methy450PP}
\title{
Data preprocessing step includes normalization, transformation,gene filtering.
%%  ~~function to do ... ~~
}
\description{
  It allows user to choose several filtering step or modify filtering criteria for specific quality control purpose. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lumi.methy450PP(data, na.omit = TRUE, normalization = FALSE, transfm = FALSE, filtersample = TRUE, filterdetectP = TRUE, detectPcut = 0.75, locidiff = FALSE, locidiffcut = 0.01, Xchrom = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
   a exprmethy450 class returned by the lumi.methy450R function
%%     ~~Describe \code{data} here~~
}
\item{na.omit}{
  if TRUE remove the sites containing na value
%%     ~~Describe \code{na.omit} here~~
}
\item{normalization}{
  if TRUE, quantile normalization performed
%%     ~~Describe \code{normalization} here~~
}
\item{transfm}{
  if TRUE, arcsin square root of beta matrix is performed
%%     ~~Describe \code{transfm} here~~
}
\item{filtersample}{
  keep the samples with detection Pvalue <10e-5 on more than 75% sites  
%%     ~~Describe \code{filtersample} here~~
}
\item{filterdetectP}{
   Remove the sites with detection Pvalue > 0.01 across  detectPcut percentage of the samples
  
%%     ~~Describe \code{filterdetectP} here~~
}
\item{detectPcut}{
Remove the sites with detection Pvalue > 0.01 across 75\% of the samples
%%     ~~Describe \code{detectPcut} here~~
}
\item{locidiff}{
  if TRUE, keep the sites has big variance between two groups.
%%     ~~Describe \code{locidiff} here~~
}
\item{locidiffcut}{
   if locidiff = TRUE, keep the sites has variance great than the locidiffcut between two groups 
%%     ~~Describe \code{locidiffcut} here~~
}
\item{Xchrom}{
  if TRUE, Remove the sites on chromosome X
%%     ~~Describe \code{Xchrom} here~~
}
}
\details{
  It allows user to choose several filtering step or modify filtering criteria for specific quality control purpose. By default, IMA will filter out loci with missing beta value, from the X chromosome or with median detection P-value greater than 0.05. The option for sample level quality control is also provided (B. C. Christensen et al., 2011). Although the raw beta values will be analyzed as recommended by Illumina, the user can choose Arcsine squareroot transformation.  Both scaling and quantile normalization options are available for cross sample normalization.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  return a methy450batch class including:
  \item{bmatrix}{ beta matrix for each ~450k site on each sample}
\item{detectP}{ a matrix of detection pvalue for each ~450k site on each sample}
\item{annot}{a matrix of annotation for each ~450k site}
\item{groupinfo}{ a list of sample ID and penotype for each sample}
  \item{TSS1500Ind}{ a list of sites ID only on the TSS1500 region for each gene  }
  \item{TSS200Ind}{ a list of sites ID only on the TSS200 region for each gene }
  \item{UTR5Ind}{ a list of sites ID only on the 5' UTR region for each gene }
  \item{EXON1Ind}{a list of sites ID only on the 1st EXON for each gene}
  \item{UTR3Ind}{a list of sites ID only on the 3' UTR region for each gene}
  \item{PROMOTERInd}{a list of sites ID only on the promoter region for each gene}
  \item{GENEBODYInd}{a list of sites ID only on the gene body region for each  gene}
  \item{GENEInd}{a list of sites ID for each gene}
  \item{ISLANDInd}{a list of istes ID only on the ISLAND region for each UCSC_CPG_ISLAND}
  \item{ NSHOREInd}{a list of istes ID only on the N Shore region for each  UCSC_CPG_ISLAND}
  \item{SSHOREInd}{a list of istes ID only on the S Shore region for each  UCSC_CPG_ISLAND}
  \item{NSHELFInd}{a list of istes ID only on the N Shelf region for each   UCSC_CPG_ISLAND}
  \item{SSHELFInd}{a list of istes ID only on the S Shelf region for each   UCSC_CPG_ISLAND}
  \item{UCSCInd}{a list of istes ID for each  UCSC_CPG_ISLAND}
  
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  https://r-forge.r-project.org/projects/ima/
%% ~put references to the literature/web site here ~
}
\author{
  Dan Wang; Li Yan; Qiang Hu; Song Liu 
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~regionswrapper }
\keyword{ ~lumi.methy450R }% __ONLY ONE__ keyword per line
