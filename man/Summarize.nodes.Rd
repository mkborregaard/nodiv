\name{Summarize.nodes}
\alias{Summarize.nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Summarize.nodes(dispersion, comm = hcom, tree = htree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dispersion}{
%%     ~~Describe \code{dispersion} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dispersion, comm = hcom, tree = htree) 
{
    nodenums <- as.numeric(colnames(dispersion))
    overdisp <- apply(dispersion, 2, function(x) mean(abs(x), 
        na.rm = T))
    overdisp[1] <- 0
    overdisp[overdisp == Inf] <- 0
    nodesize <- sapply(nodenums, Node.size)
    nodeoccupancy <- sapply(nodenums, Node.occupancy, comm = comm, 
        tree = tree)
    nodestats <- data.frame(nodenums, overdisp, nodeoccupancy, 
        nodesize)
    return(nodestats)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
