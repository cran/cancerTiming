\name{eventTiming}
\alias{eventTiming}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the time of events in tumor data
}
\description{
Estimate the proportion of time spent between different chromosomal abnormalities based on the allele frequencies of mutated locations. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eventTiming(x, m, history, totalCopy, method=c("fullMLE","partialMLE","Bayes"),type=c("gain","CNLOH"),seqError=0, bootstrapCI=NULL, B=if(method=="Bayes") 10000 else 500, CILevel=0.95, alpha=1,tdf=4,normCont=0,verbose=TRUE,returnAssignments=FALSE,coverageCutoff=1,minMutations=10,returnData=FALSE,init=NULL,maxiter=100, tol=0.0001)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector. the number of reads/fragments containing the variant
%%     ~~Describe \code{x} here~~
}
  \item{m}{ vector. the number of reads/fragments covering the location with the variant (the coverage)
%%     ~~Describe \code{m} here~~
}
  \item{history}{ a matrix, based on the history of the region (see Details)
%%     ~~Describe \code{history} here~~
}
  \item{totalCopy}{ integer. the total number of copies of the tumor DNA for this region 
%%     ~~Describe \code{totalCopy} here~~
}
\item{method}{what estimation method to use, one of ``fullMLE",``partialMLE",``Bayes" }
  \item{type}{ type of region, either a gain or a CNLOH region
%%     ~~Describe \code{type} here~~
}
\item{seqError}{Probability of sequencing error}
  \item{bootstrapCI}{ type of bootstrap confidence interval to calculate, one of ``parametric", ``nonparametric". If NULL, then the confidence interval is not calculated
%%     ~~Describe \code{bootstrapCI} here~~
}
  \item{B}{ number of bootstrap samples to take
%%     ~~Describe \code{B} here~~
}
\item{CILevel}{
	At what level the confidence intervals should be calculated.
}
 \item{alpha}{parameter the Dirichlet prior of the bayesian estimate}
\item{tdf}{parameter for the number of degrees of freedom for the t proposal density used in the bayesian estimate}
  \item{normCont}{ the proportion of normal contamination, between 0 and 1.
%%     ~~Describe \code{normCont} here~~
}
  \item{verbose}{ logical. Whether to give additional warnings as the program is running. 
%%     ~~Describe \code{verbose} here~~
}
  \item{returnAssignments}{ logical. Whether to re
%%     ~~Describe \code{returnAssignments} here~~
}
  \item{coverageCutoff}{ minimum value for m[i]; any entries with m[i]<coverageCutoff are ignored in estimation.
%%     ~~Describe \code{coverageCutoff} here~~
}
  \item{minMutations}{ minimum number of mutations required.
%%     ~~Describe \code{minMutations} here~~
}
  \item{returnData}{ logical. Whether to return the input data, x and m, as a data.frame
%%     ~~Describe \code{returnData} here~~
}
  \item{init}{ initial value of multinomial parameter q passed to estimateQ. 
%%     ~~Describe \code{init} here~~
}
  \item{maxiter}{maximum number of iterations in calculation q.
%%     ~~Describe \code{maxiter} here~~
}
  \item{tol}{tolerance in the convergence of q
%%     ~~Describe \code{tol} here~~
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
A list with values
 \item{pi }{estimate of pi vector}
 \item{piCI }{bootstrap confidence interval, if requested}
\item{q}{estimate of the multinomial parameter q}
\item{perLocationProb}{output from estimateQ giving per location P(P[i]|X[i],q), if requested. Only locations used in the estimation are included.}
\item{optimDetails}{optimization details from estimateQ}
\item{call}{list of the parameters of the call to the function: history,totalCopy,type,exactAllele,normCont,coverageCutoff,minMutations. In addition, `alleleSet' is included in this list, which is the set of possible alleles for this history, after adjusting for normal contamination. }
}
\references{
Durinck, S, et al. (2011). ``Temporal Dissection of Tumorigenesis in Primary Cancers." Cancer Discovery, 1(2): 137-143.
Greenman, C D et al. (2012). ``Estimation of rearrangement phylogeny for cancer genomes." Genome Research, 22(2): 346-361.
Purdom, E, et al. (submitted). ``Timing Chromosomal Abnormalities within Cancer Samples.''
}
\author{
Elizabeth Purdom
}

\examples{
data(mutData)
ACNLOH<-matrix(c(1,3,1,0),ncol=2,nrow=2,byrow=TRUE)
onlyMuts<-subset(mutData,is.na(rsID) & position <= 1.8E7)
onlyMuts$t_depth<-onlyMuts$t_ref_count+onlyMuts$t_alt_count
x<-eventTiming(x=onlyMuts$t_alt_count,m=onlyMuts$t_depth,history=ACNLOH,totalCopy=2,type="CNLOH",normCont=0.22)
x$pi #estimate of time of stages
x$q #estimate of the multinomial (likelihood of each of the alleles)
x$call$alleleSet #possible set of alleles after adjusting for normal contamination
}
