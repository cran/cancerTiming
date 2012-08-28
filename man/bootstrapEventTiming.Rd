\name{bootstrapEventTiming}
\alias{bootstrapEventTiming}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap the results of eventTiming
}
\description{
Create bootstrap estimates of pi base on the results of a call to eventTiming.
}
\usage{
bootstrapEventTiming(B,  type = c("parametric", "nonparametric"),call, pi, x, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{ number of bootstrap samples to take
%%     ~~Describe \code{B} here~~
}
  \item{type}{type of bootstrap confidence interval to calculate, one of ``parametric", ``nonparametric"
%%     ~~Describe \code{type} here~~
}
  \item{call}{
list with elements (given as `call' in the output of eventTiming): ``alleleSet", ``history", ``totalCopy", ``type", ``exactAllele", ``normCont", ``coverageCutoff", ``minMutations", ``init", ``maxiter", ``tol"
}
  \item{pi}{
the estimate of $pi$ from which to create bootstrap samples (if type="parametric" )
}
  \item{x}{ vector. the number of reads/fragments containing the variant; not needed if ``inputData'' included as part of call
%%     ~~Describe \code{x} here~~
}
  \item{m}{ vector. the number of reads/fragments covering the location with the variant (the coverage)
%%     ~~Describe \code{m} here~~
}

}
\value{
Matrix of dimension (B,length(pi)) with the estimate of pi for each of the bootstrap samples
}
\author{
Elizabeth Purdom}
\examples{
##can do this within eventTiming function, but here is an example doing it separately...
data(mutData)
ACNLOH<-matrix(c(1,3,1,0),ncol=2,nrow=2,byrow=TRUE)
onlyMuts<-subset(mutData,is.na(rsID) & position <= 1.8E7)
onlyMuts$t_depth<-onlyMuts$t_ref_count+onlyMuts$t_alt_count
x<-eventTiming(x=onlyMuts$t_alt_count,m=onlyMuts$t_depth,history=ACNLOH,totalCopy=2,type="CNLOH",normCont=0.22)
piBoot<-bootstrapEventTiming(B=100,pi=x$pi,call=x$call,x=onlyMuts$t_alt_count,m=onlyMuts$t_depth,type="parametric")
}
