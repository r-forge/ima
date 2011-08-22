\name{sigtest}
\alias{sigtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Differential Methylation Analysis
%%  ~~function to do ... ~~
}
\description{
  For each specific region, Wilcoxon rank-sum test (default), Student's t-test and empirical Bayes statistics (G. K. Smyth, 2004) are available for inference in differential testing. Robust linear regression is available as an option to infer methylation change associated with continuous variable (e.g., age). A variety of multiple testing correction algorithms is available, including stringent Bonferroni correction and widely used false discovery rate control. Users can specify the significance criteria in parameter file. The same statistical inference and multiple test correction procedures described above will also be applied to each single site to obtain site-level differential methylation inference.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigtest(dataf, gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = "BH", concov = "OFF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a  methy450batch class returned by the lumi.methy450PP function
    
%%     ~~Describe \code{dataf} here~~
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
