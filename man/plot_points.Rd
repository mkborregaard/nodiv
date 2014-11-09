\name{plot_points}
\alias{plot_points}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_points(x, coords, pch = 16, col = rev(terrain.colors(255)), shape = NULL, shapefill = "grey", shapeborder = "grey", zlims, bg = par("bg"), zoom_to_points = FALSE, xlim = NULL, ylim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{shape}{
%%     ~~Describe \code{shape} here~~
}
  \item{shapefill}{
%%     ~~Describe \code{shapefill} here~~
}
  \item{shapeborder}{
%%     ~~Describe \code{shapeborder} here~~
}
  \item{zlims}{
%%     ~~Describe \code{zlims} here~~
}
  \item{bg}{
%%     ~~Describe \code{bg} here~~
}
  \item{zoom_to_points}{
%%     ~~Describe \code{zoom_to_points} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
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
function (x, coords, pch = 16, col = rev(terrain.colors(255)), 
    shape = NULL, shapefill = "grey", shapeborder = "grey", zlims, 
    bg = par("bg"), zoom_to_points = FALSE, xlim = NULL, ylim = NULL, 
    ...) 
{
    if (inherits(x, "SpatialPointsDataFrame")) {
        x <- x@data
        coords <- SpatialPoints(x)
    }
    else {
        if (missing(coords)) 
            stop("coords must be defined if x is not a SpatialPixelsDataFrame")
        else {
            if (!inherits(coords, "SpatialPoints")) {
                if (is.matrix(coords)) 
                  coords <- as.data.frame(coords)
                if (is.data.frame(coords)) {
                  if (!is.vector(x)) 
                    stop("x must be a vector or a SpatialPointsDataFrame")
                  if (!nrow(coords) == length(x)) 
                    stop("x and coords must have the same number of elements")
                  if (!ncol(coords) == 2) 
                    stop("coords must have exactly two columns, for the x and y coordinates of data")
                  coords <- SpatialPoints(coords)
                }
                else stop("Undefined arguments")
            }
        }
    }
    if (zoom_to_points & !is.null(shape)) {
        if (is.null(xlim)) 
            xlim <- 1.1 * bbox(coords)[1, ]
        if (is.null(ylim)) 
            ylim <- 1.1 * bbox(coords)[2, ]
    }
    if (missing(zlims)) 
        zlims <- c(min(x, na.rm = T), max(x, na.rm = T))
    coords <- SpatialPoints(coords)
    oldpar <- par()
    par(mar = c(5, 4, 4, 6) + 0.1)
    if (is.null(shape)) 
        plot(coords, col = create.cols(x, col, zlim = zlims), 
            pch = pch, ...)
    else {
        plot(shape, col = shapefill, border = shapeborder, xlim = xlim, 
            ylim = ylim, ...)
        plot(coords, col = create.cols(x, col, zlim = zlims), 
            pch = pch, bg = bg, add = T)
    }
    library(fields)
    par <- oldpar
    image.plot(zlim = zlims, legend.only = TRUE, smallplot = c(0.85, 
        0.87, 0.38, 0.65), col = col)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
