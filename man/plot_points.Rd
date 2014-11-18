\name{plot_points}
\alias{plot_points}
\alias{plot_grid}
\alias{plot_sitestat}

\title{
Plot a variable on spatial points or a spatial grid using a color scale.
}
\description{
Plots a variable either as colored points or as a colored raster. \code{plot_sitestat} is a wrapper function that dispatches the other functions according to the \code{type} of \code{distrib_data}.
}
\usage{
plot_points(x, coords, col = rev(terrain.colors(64)), shape = NULL, 
  shapefill = "grey",  zlim= NULL,  zoom_to_points = FALSE, 
  pch = 16, bg = par("bg"), ...)
plot_grid(x, coords, col = rev(terrain.colors(64)), shape = NULL, 
  shapefill = "grey",  zlim = NULL, 
  zoom_to_points = FALSE,  ...)
plot_sitestat(distrib_data, x, ...)
}

\arguments{
  \item{x}{
A numerical variable to be plotted using a color scale
}
  \item{coords}{
A matrix or \code{data.frame} with two columns, or an object inheriting from \code{SpatialPoints}, of the same length as \code{x}
}
  \item{col}{
The color palette to use for the color scale
}
  \item{shape}{
A shape file giving a map to use as a background for the plot
}
  \item{shapefill}{
If shape is a polygon, the color to use to fill the polygons. If shape is a line shapefile, the color to use for plotting lines
}
  \item{zlim}{
The values to use for the extremes of the color scale
}
  \item{zoom_to_points}{
Given a shapefile, the function will set the outer limits of the plot to the bounding box of the shapefile. If \code{zoom_to_points} is \code{TRUE}, the outer limits will be set by the bounding box of the data points instead.
}
  \item{pch}{
The point character to use for plotting points
}
  \item{bg}{
The central color of points when using a pch value between 21 and 25
}
  \item{distrib_data}{
An object of type \code{distrib_data}, \code{nodiv_data} or \code{nodiv_result} to use for plotting x.
}
  \item{\dots}{
Additional arguments to pass to plot
}
}

\author{
Michael Krabbe Borregaard 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
