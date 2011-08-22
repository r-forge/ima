\name{sumregionfun}
\alias{sumregionfun}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Methylation Index Calculation
%%  ~~function to do ... ~~
}
\description{
  for each specific region (e.g., promoter), IMA will collect the loci within it and derive an index of overall region methylation value. Currently, there are three different index metrics implemented in IMA: mean, median, and Tukey's Biweight robust average.  By default, the median beta values will be used as the region's methylation index for further analysis. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sumregionfun(indexlist = dataf@GENEInd, beta, sumregion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexlist}{specify which annotated region you are interested
%%     ~~Describe \code{indexlist} here~~
}
\item{beta}{
  beta matrix for each ~450k site
%%     ~~Describe \code{beta} here~~
}
\item{sumregion}{
  The method derive the overall region methylation value for each gene or UCSC_CPG_ISLAND. The options are mean, median, and Tukey's Biweight robust average.
%%     ~~Describe \code{sumregion} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  return a summarized matrix for the annotated region
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
