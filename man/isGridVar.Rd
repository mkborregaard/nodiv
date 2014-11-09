\name{isGridVar}
\alias{isGridVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
isGridVar(gridVar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gridVar}{
%%     ~~Describe \code{gridVar} here~~
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
function (gridVar) 
{
    dists <- diff(sort(unique(gridVar)))
    distab <- table(dists)
    smallest <- as.numeric(names(distab[1]))
    most_common <- as.numeric(names(distab))[distab == max(distab)]
    return(isTRUE(all.equal(dists/smallest, floor(dists/smallest))) & 
        smallest \%in\% most_common)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
