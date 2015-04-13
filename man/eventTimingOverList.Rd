\name{eventTimingOverList}
\alias{eventTimingOverList}
\alias{getPi0Summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
eventTiming for multiple samples and regions
}
\description{
eventTimingOverList is a wrapper to eventTiming that allows for timing of common events over several regions of a sample and/or multiple samples. getPi0Summary gets summary information about \eqn{\pi_0}{pi[0]} (the first event) from the output of that function and returns a simple dataframe of the estimate of \eqn{\pi_0}{pi[0]} for every region and sample combination. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eventTimingOverList(dfList, normCont, eventArgs)
getPi0Summary(eventList, CI = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfList}{a list of mutation data per sample. Each element of the list corresponds to a sample, and consists of a data frame of mutation information for all of the regions that are to be timed. The data frame must have a certain format, see description below.
}
  \item{normCont}{
a vector of equal length as dfList giving the normal contamination of each sample
}
  \item{eventArgs}{
	list of arguments passed to eventTiming via `do.call'. Should NOT contain the arguments `x',`m',`history',`totalCopy',`type',`mutationId' or `normCont'
%%     ~~Describe \code{eventArgs} here~~
}
  \item{eventList}{
Output of eventTimingOverList (see below)
}
  \item{CI}{logical, whether to grab CI from the output. 
%%     ~~Describe \code{CI} here~~
}
}
\details{
The data frame of mutation data for each sample must have the following columns: `segId', `type', `nMutAllele', `nReads', `mutationId'. The rows of the data frame correspond to mutated locations within the sample. `segId' is an identifier of the segmented region that the mutation is in; eventTiming will be run using the mutation data for each `segId' in the sample. `type' identifies the type of segment, one of c("Other","CNLOH","SingleGain","Diploid","DoubleGain"). `nMutAllele' gives for every row (mutated location) the number of reads with the mutation and corresponds to `x' imputed into eventTiming. `nReads' gives the total number of reads covering the location and corresponds to `m' imputed into eventTiming. `mutationId' is a unique identifier for the mutation. 

For arguments passed to eventTiming, it is recommended to pass verbose=FALSE if you want to avoid many warnings about mutations that do not meet the necessary criteria of coverage, etc.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
eventTimingOverList returns a list equal to the number of samples. Each element of that list (i.e. per sample) is itself a list of length three corresponding to the three types of regions that can be timed, ``CNLOH",``SingleGain", and ``DoubleGain". Each of these gives the output of eventTiming per region of that type.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Elizabeth Purdom
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{eventTiming}}
}
