\name{plot_richness}
\alias{plot_richness}

\title{
Plot a richness map for an object of type distrib_data, nodiv_data or nodiv_result
}
\description{
Create a raster or point plot showing the species richness of an object. Useful for comparing the species richness of sister clades
}
\usage{
plot_richness(distrib_data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distrib_data}{
An object of types \code{distrib_data}, \code{nodiv_data} or \code{nodiv_result}
}
  \item{\dots}{
Further arguments to pass to the plotting function
}
}
\details{
The plot will be done as a spatial grid or as colored points, depending on the \code{type} of the \code{nodiv_result} object
}
\value{
If \code{distrib_data} has \code{type raster}, a raster object of the SOS values is returned. 
}

\author{
Michael Krabbe Borregaard
}


%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
