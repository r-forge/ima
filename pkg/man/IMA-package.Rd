\name{IMA-package}
\alias{IMA-package}

\docType{package}
\title{
IMA(Illumina Methylation Analyzer) is a package designed to automate the pipeline for analyzing site-level and region-level methylation changes in epigenetic studies utilizing the 450K DNA methylation microarray}
\description{
IMA automates the tasks commonly required for the differential analysis of epigenetic data sets utilizing the 450K DNA methylation microarray. The package makes use of Illumina methylation annotation for region definition, as well as several Bioconductor packages for various preprocessing and differential testing steps. Written in open source R environment, it provides the flexibility for users to adopt, extend and customize the functionality for their specific needs. It can be used as an automatic pipeline to analyze specific regions as well as specific sites for downstream functional exploration and hypothesis generation.
}

\details{
\tabular{ll}{
Package: \tab ima\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-25\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Dan Wang, Li Yang, Qiang Hu,Dominic J Smiraglia,Song Liu
Maintainer: Dan Wang <wangdan412@gmail.com>
}
\seealso{
  \code{\link{IMA.methy450R}},\code{\link{IMA.methy450PP}},\code{\link{regionswrapper}},\code{\link{testfunc}},\code{\link{sigtest}},\code{\link{sumregionfun}}
}
\examples{
###Users specify the data paths/paramters in the pipeline.R before
#run the comman below 
###R --no-save < pipeline.R##This will auto-mate to produce the result
#of methylation
###the pipeline.R file could be found here
###http://ima.r-forge.r-project.org/
}
