\name{SOS}
\alias{SOS}
\alias{plotSOS}

\title{
Gives the SOS values for a given node.
}
\description{
\code{SOS} gives the SOS value for the node, used in the interpretation of node-based analysis of species distributions after running \code{Node_analysis}.  \code{plotSOS} plots the values on a map
}
\usage{
SOS(nodiv_result, node)
plotSOS(nodiv_result, node, col = cm.colors(64), zlim, ...)
}

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
The plot will be done as a spatial grid or as colored points, depending on the \code{type} of the \code{nodiv_result} object. If the \code{nodiv_result} object has a \code{shape} object of type \code{raster}, this will be plotted normally in the background - but resizing the window will cause the plots/pixels and the underlying raster object to lose alignment, due to a bug in the \code{raster} package.
}
\value{
\code{SOS} returns a numerical vector with the SOS values. If \code{nodiv_result} has \code{type raster}, a raster object of the SOS values is returned by \code{plotSOS}. 
}
\references{
Borregaard, M.K., Rahbek, C., Fjeldsaa, J., Parra, J.L., Whittaker, R.J. & Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12283
}
\author{
Michael Krabbe Borregaard
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Node_analysis}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
