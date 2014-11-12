\name{plotSOS}
\alias{plotSOS}

\title{
Plot a map of SOS values for a given node
}
\description{
Create a plot used in the interpretation of node-based analysis of species distributions, after running \code{Node_analysis}
}
\usage{
plotSOS(nodiv_result, node, col = cm.colors(64), zlim, ...)
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
