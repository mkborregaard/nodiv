\name{rdtable_nodesig}
\alias{rdtable_nodesig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rdtable_nodesig(simcom, node.sp, repeats, show = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simcom}{
%%     ~~Describe \code{simcom} here~~
}
  \item{node.sp}{
%%     ~~Describe \code{node.sp} here~~
}
  \item{repeats}{
%%     ~~Describe \code{repeats} here~~
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
function (simcom, node.sp, repeats, show = F) 
{
    tempcom <- cbind(rowSums(simcom[, node.sp]), rowSums(simcom[, 
        !node.sp]))
    r <- rowSums(tempcom)
    c <- colSums(tempcom)
    ll <- 0
    nodereps <- replicate(repeats - 1, {
        simcom <- r2dtable(1, r, c)[[1]]
        if (show) 
            print(ll <<- ll + 1)
        simcom[, 1]
    })
    return(nodereps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
