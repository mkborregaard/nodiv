\name{Nodesig}
\alias{Nodesig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Nodesig(commat, node.sp, repeats = 100, method = c("quasiswap", "rdtable"), show = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{commat}{
%%     ~~Describe \code{commat} here~~
}
  \item{node.sp}{
%%     ~~Describe \code{node.sp} here~~
}
  \item{repeats}{
%%     ~~Describe \code{repeats} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{show}{
%%     ~~Describe \code{show} here~~
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
function (commat, node.sp, repeats = 100, method = c("quasiswap", 
    "rdtable"), show = F) 
{
    if (sum(node.sp) == 1 | sum(!node.sp) == 1) 
        return(NA)
    method = match.arg(method)
    require(vegan)
    require(fields)
    simcom <- commat
    nodereps <- switch(method, quasiswap = quasiswap_nodesig(simcom, 
        node.sp, repeats, show), rdtable = rdtable_nodesig(simcom, 
        node.sp, repeats, show))
    nodeemp <- rowSums(commat[, node.sp])
    nodereps <- cbind(nodeemp, nodereps)
    ord <- apply(nodereps, 1, rank)[1, ]
    ord[ord == repeats] <- repeats - 1
    ord[rowSums(commat) == 0] <- NA
    nodemeans <- rowMeans(nodereps)
    nodesds <- apply(nodereps, 1, sd)
    ses <- (nodeemp - nodemeans)/nodesds
    pval <- ord/repeats
    return(list(ses = ses, pval = pval, nodeemp = nodeemp, nodemeans = nodemeans, 
        nodesds = nodesds))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
