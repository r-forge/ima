\name{testfunc}
\alias{testfunc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Differential testing function
%%  ~~function to do ... ~~
}
\description{
  pick up the test for differential methylation testing 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
testfunc(eset = eset, concov = c("ON", "OFF"), sumregion = c("mean", "median", "tbrm"), testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), lev1 = lev1, lev2 = lev2, grouplev = grouplev)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
    a beta matrix
%%     ~~Describe \code{eset} here~~
}
\item{concov}{
  if "ON", covariates is continuous variable
%%     ~~Describe \code{concov} here~~
}
\item{testmethod}{
   Wilcoxon rank-sum test (default), Student's t-test and empirical Bayes statistics (G. K. Smyth, 2004).
%%     ~~Describe \code{testmethod} here~~
}
\item{Padj}{
  specify which method applied for multiple testing correction.
%%     ~~Describe \code{Padj} here~~
}
\item{lev1}{
  specify which phenotype group to be used for comparison
%%     ~~Describe \code{lev1} here~~
}
\item{lev2}{
  specify which phenotype group to be used for comparison
%%     ~~Describe \code{lev2} here~~
}
\item{grouplev}{
  specify the two phenotype groups to be used for comparison
%%     ~~Describe \code{grouplev} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  return a matrix with "P-Value", "Adjust Pval","beta-Difference" in columns for each sites or annotated region
  
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
