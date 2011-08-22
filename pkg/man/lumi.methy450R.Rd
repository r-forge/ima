\name{lumi.methy450R}
\alias{lumi.methy450R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Readin methy450k data
%%  ~~function to do ... ~~
}
\description{
This function is used to read the methy450k data into R. This require two input file, one is the txt file containing all the methylation information and produced from the illumina BeadStudio or GenomeStudio
software. The other one is the sample.txt file which specify the penotype for each samples 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lumi.methy450R(file = "../data/SampleMethFinalReport_Smiraglia.txt", columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"), groupfile = "../data/sample.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
 this file contains all the methylation
information and was produced from the illumina BeadStudio or GenomeStudio
software.
%%     ~~Describe \code{file} here~~
}
\item{columnGrepPattern}{
 specify the string pattern's name of beta value and detection pvalue in the .txt file 
%%     ~~Describe \code{columnGrepPattern} here~~
}
\item{groupfile}{
 This  file specify the penotype for each samples 
%%     ~~Describe \code{groupfile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It will return a exprmethy450 class including below values:a beta matrix, a annot matrix, a matrix with detect_pvalue for each ~450k site and a group info
matrix.
\item{bmatrix}{ beta matrix for each ~450k site on each sample}
\item{detectP}{ a matrix of detection pvalue for each ~450k site on each sample}
\item{annot}{a matrix of annotation for each ~450k site}
\item{groupinfo}{ a list of sample ID and penotype for each sample} 
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
\keyword{~lumi.methy450PP}
\keyword{~regionswrapper}
