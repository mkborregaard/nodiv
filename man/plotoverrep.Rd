\name{plotoverrep}
\alias{plotoverrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotoverrep(node_number, dispersion_eff = rep_matrix, coords = dat.LL, tree = htree, plottype = c("map", "points"), new.window = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_number}{
%%     ~~Describe \code{node_number} here~~
}
  \item{dispersion_eff}{
%%     ~~Describe \code{dispersion_eff} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{plottype}{
%%     ~~Describe \code{plottype} here~~
}
  \item{new.window}{
%%     ~~Describe \code{new.window} here~~
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
function (node_number, dispersion_eff = rep_matrix, coords = dat.LL, 
    tree = htree, plottype = c("map", "points"), new.window = FALSE, 
    ...) 
{
    plottype = match.arg(plottype)
    if (plottype == "map") {
        overrep <- dispersion_eff[, which(colnames(dispersion_eff) == 
            as.character(node_number))]
        max_z = max(abs(overrep[abs(overrep) != Inf]), na.rm = T)
        zlims = c(-max_z, max_z)
        map.var(overrep, coords$Long, coords$Lat, new.window, 
            zlim = zlims, ...)
    }
    else {
        overrep <- dispersion_eff[, which(nodenumbers(tree) == 
            node_number)]
        oldpar <- par()
        par(mar = c(5, 4, 4, 4) + 0.1)
        max_z = max(abs(overrep), na.rm = T)
        zlims = c(-max_z, max_z)
        cols = create.cols(overrep, zlims = zlims)
        plot(coords$Long, coords$Lat, col = cols, pch = 16, cex = 1.5, 
            ...)
        plot(sa, add = T)
        library(fields)
        par <- oldpar
        image.plot(zlim = zlims, col = rev(rainbow(256, start = 0.02, 
            end = 0.59)), legend.only = T)
    }
    title(paste("representation eff size of node", node_number))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
