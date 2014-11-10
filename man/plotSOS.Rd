\name{plotSOS}
\alias{plotSOS}

\title{
Plot a map of SOS values for a given node, after running \code{Node_analysis}
}
\description{
Create a plot used in the interpretation of node-based analysis of species distributions
}
\usage{
plotSOS(nodiv_result, node, col = brewer.pal(11, "RdYlBu"), zlim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodiv_result}{
The result of running \code{Node_analysis}
}
  \item{node}{
The node label or node number specifying the node to plot
}
  \item{col}{
The color palette to use for illustrating SOS values
}
  \item{zlim}{
The values to use for the extremes of the color palette
}
  \item{\dots}{
Further arguments to pass to the plotting function
}
}
\details{
The plot will be done as a spatial grid or as colored points, depending on the \code{type} of the \code{nodiv_result} object
}
\value{
If \code{nodiv_result} has \code{type raster}, a raster object of the SOS values is returned. 
}
\references{
Borregaard, M.K., Rahbek, C., Fjeldså, J., Parra, J.L., Whittaker, R.J. & Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12283
}
\author{
Michael Krabbe Borregaard
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Node_analysis}}
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
        plot_grid(SOS, nodiv_result$coords, col = brewer.pal(11, 
            "RdYlBu"), zlim = zlim, shape = shape, ...)
    else plot_points(SOS, nodiv_result$coords, col = brewer.pal(11, 
        "RdYlBu"), zlim = zlim, shape = shape, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
