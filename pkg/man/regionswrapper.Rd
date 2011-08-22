\name{regionswrapper}
\alias{regionswrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test on the 14 annotated region  
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regionswrapper(dataf, sumregion = c("mean", "median", "tbrm"), gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = "BH", concov = c("OFF", "ON"), list14excel = list14excel, list14Rdata = list14Rdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
     a  methy450batch class returned by the lumi.methy450PP function
    
%%     ~~Describe \code{dataf} here~~
}
\item{sumregion}{
  
%%     ~~Describe \code{sumregion} here~~
}
\item{gcase}{
   the string name of case group in the sample.txt file
%%     ~~Describe \code{gcase} here~~
}
\item{gcontrol}{
  the string name of control group in the sample.txt file
%%     ~~Describe \code{gcontrol} here~~
}
\item{testmethod}{
  Wilcoxon rank-sum test (default), Student's t-test and empirical Bayes statistics (G. K. Smyth, 2004).
%%     ~~Describe \code{testmethod} here~~
}
\item{Padj}{
   specify which method applied for multiple testing correction.
%%     ~~Describe \code{Padj} here~~
}
\item{concov}{
   if "ON", covariates is continuous variable
%%     ~~Describe \code{concov} here~~
}
\item{list14excel}{
the name of output excel file containing the differential testing result for each annotated region
%%     ~~Describe \code{list14excel} here~~
}
\item{list14Rdata}{
the name of Rdata containing the differential testing result for each annotated region
%%     ~~Describe \code{list14Rdata} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return an excel file containing the differential testing result for each annotated region and a Rdata containing the differential testing  result for each annotated region
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
