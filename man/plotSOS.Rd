\name{plotSOS}
\alias{plotSOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSOS(nodiv_result, node, col = brewer.pal(11, "RdYlBu"), zlim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodiv_result}{
%%     ~~Describe \code{nodiv_result} here~~
}
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{zlim}{
%%     ~~Describe \code{zlim} here~~
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
function (nodiv_result, node, col = brewer.pal(11, "RdYlBu"), 
    zlim, ...) 
{
    if (!inherits(nodiv_result, "nodiv_result")) 
        stop("nodiv_result must be the result object from running Node_analysis")
    node <- node[1]
    if (is.character(node)) {
        if (is.null(nodiv_result$phylo$node.label)) 
            stop("node could not be matched, as the phylogeny does not have node labels")
        node <- match(node, nodiv_result$phylo$node.label)
        if (is.na(node)) 
            stop("the node could not be matched to the node labels")
    }
    if (node > Nspecies(nodiv_result)) 
        node <- node - Nspecies(nodiv_result)
    SOS <- nodiv_result$SOS[, node]
    if (missing(zlim)) {
        maxabs <- max(abs(SOS))
        zlim <- c(-maxabs, maxabs)
    }
    if (is.null(nodiv_result$shape)) 
        shape <- NULL
    else shape <- nodiv_result$shape
    if (nodiv_result$type == "grid") 
        map_var(SOS, nodiv_result$coords, col = brewer.pal(11, 
            "RdYlBu"), zlim = zlim, shape = shape, ...)
    else plot_points(SOS, nodiv_result$coords, col = brewer.pal(11, 
        "RdYlBu"), zlim = zlim, shape = shape, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
