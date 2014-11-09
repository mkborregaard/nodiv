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
Nodesig(nodiv_data, Node_sp, repeats = 100, method = c("quasiswap", "rdtable"), show = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodiv_data}{
%%     ~~Describe \code{nodiv_data} here~~
}
  \item{Node_sp}{
%%     ~~Describe \code{Node_sp} here~~
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
function (nodiv_data, Node_sp, repeats = 100, method = c("quasiswap", 
    "rdtable"), show = F) 
{
    if (sum(Node_sp) == 1 | sum(!Node_sp) == 1) 
        return(rep(NA, 5))
    method = match.arg(method)
    require(vegan)
    require(fields)
    simcom <- nodiv_data$comm
    nodereps <- switch(method, quasiswap = quasiswap_nodesig(simcom, 
        Node_sp, repeats, show), rdtable = rdtable_nodesig(simcom, 
        Node_sp, repeats, show))
    nodeemp <- rowSums(nodiv_data$comm[, Node_sp])
    nodereps <- cbind(nodeemp, nodereps)
    ord <- apply(nodereps, 1, rank)[1, ]
    ord[ord == repeats] <- repeats - 1
    ord[rowSums(nodiv_data$comm) == 0] <- NA
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
