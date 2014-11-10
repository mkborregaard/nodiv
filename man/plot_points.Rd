\name{plot_points}
\alias{plot_points}
\alias{plot_grid}

\title{
Plot a variable on spatial points or a spatial grid using a color scale
}
\description{
Plots the variable either as colored points or as a colored raster
}
\usage{
plot_points(x, coords, col = rev(terrain.colors(255)), shape = NULL, shapefill = "grey", shapeborder = "grey", zlim= c(min(vec,na.rm = T),max(vec,na.rm = T)),  zoom_to_points = FALSE, xlim = NULL, ylim = NULL, pch = 16, bg = par("bg"), ...)
plot_grid(x, coords, col = rev(terrain.colors(255)), shape = NULL, shapefill = "grey", shapeborder = "grey", zlim = c(min(vec,na.rm = T),max(vec,na.rm = T)), zoom_to_points = FALSE, xlim = NULL, ylim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{shapeborder}{
If shape is a polygon, the color of the border around polygons
}
  \item{zlims}{
The values to use for the extremes of the color scale
}
  \item{zoom_to_points}{
Given a shapefile, the function will set the outer limits of the plot to the bounding box of the shapefile. If \code{zoom_to_points} is \code{TRUE}, the outer limits will be set by the bounding box of the data points instead.
}
  \item{xlim}{
The x limits of the plotting area
}
  \item{ylim}{
The y limits of the plotting area
}
  \item{pch}{
The point character to use for plotting points
}
  \item{bg}{
The central color of points when using a pch value between 21 and 25
}
  \item{\dots}{
Additional arguments to pass to plot
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{plot_grid} invisibly returns the \code{raster} object used for plotting 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Krabbe Borregaard 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
