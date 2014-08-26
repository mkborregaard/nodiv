\name{plotnode_maps_new}
\alias{plotnode_maps_new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotnode_maps_new(node_number, dispersion_eff = rep_matrix, tree = htree, comm = hcom, coords = dat.LL, new.window = FALSE, plottype = c("map", "points"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_number}{
%%     ~~Describe \code{node_number} here~~
}
  \item{dispersion_eff}{
%%     ~~Describe \code{dispersion_eff} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{new.window}{
%%     ~~Describe \code{new.window} here~~
}
  \item{plottype}{
%%     ~~Describe \code{plottype} here~~
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
function (node_number, dispersion_eff = rep_matrix, tree = htree, 
    comm = hcom, coords = dat.LL, new.window = FALSE, plottype = c("map", 
        "points")) 
{
    require(ape)
    if (new.window) 
        quartz()
    par(mfrow = c(2, 2))
    plotoverrep(node_number, dispersion_eff, coords, tree, new.window = FALSE, 
        plottype = plottype)
    plotnode(node_number, comm, tree, coords, new.window = FALSE, 
        plottype = plottype, main = paste("richness of node", 
            node_number))
    plotnode(Parent(node_number), comm, tree, coords, new.window = FALSE, 
        plottype = plottype, main = paste("richness of parent node", 
            Parent(node_number)))
    plotnode(basal.node(htree), comm, tree, coords, new.window = FALSE, 
        plottype = plottype, main = "richness of basal node")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
