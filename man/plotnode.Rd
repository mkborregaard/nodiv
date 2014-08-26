\name{plotnode}
\alias{plotnode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotnode(node_number, comm, tree = htree, coords = dat.LL, plottype = c("map", "points"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_number}{
%%     ~~Describe \code{node_number} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{plottype}{
%%     ~~Describe \code{plottype} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (node_number, comm, tree = htree, coords = dat.LL, plottype = c("map", 
    "points"), ...) 
{
    plottype = match.arg(plottype)
    if (plottype == "map") {
        richs <- Node.richness(node_number, comm, tree, coords)
        map.var(richs$richness, richs$Long, richs$Lat, ...)
    }
    else {
        richs <- Node.richness(node_number, comm, tree, coords)
        oldpar <- par()
        par(mar = c(5, 4, 4, 4) + 0.1)
        plot(richs$Long, richs$Lat, col = create.cols(richs$richness), 
            pch = 16, cex = 1.5, ...)
        plot(sa, add = T)
        library(fields)
        par <- oldpar
        image.plot(zlim = range(richs$richness, na.rm = T), col = rev(rainbow(256, 
            start = 0.02, end = 0.59)), legend.only = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
