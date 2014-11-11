\name{add_shape}
\alias{add_shape}
\title{
Add a shapefile to a \code{nodiv} data object 
}
\description{
Takes a shapefile, e.g. as points or polygons, and adds it to the data object
}
\usage{
add_shape(distrib_data, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distrib_data}{
An object of type \code{distrib_data}, \code{nodiv_data} or \code{nodiv_result}
}
  \item{shape}{
A GIS shapefile, from e.g. packages \code{maptools} or \code{rgdal}
}
}

\value{
Returns the object with the new shapefile
}

\author{
Michael Krabbe Borregaard
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{distrib_data}}
}
\examples{
## usage: data <- add_shape(data, shapefile)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
