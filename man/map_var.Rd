\name{map_var}
\alias{map_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
map_var(x, coords, shape = NULL, shapefill = "grey", shapeborder = "grey", zlim = c(min(vec, na.rm = T), max(vec, na.rm = T)), col = rev(terrain.colors(255)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
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
  \item{zlim}{
%%     ~~Describe \code{zlim} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
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
function (x, coords, shape = NULL, shapefill = "grey", shapeborder = "grey", 
    zlim = c(min(vec, na.rm = T), max(vec, na.rm = T)), col = rev(terrain.colors(255)), 
    ...) 
{
    if (inherits(x, "SpatialPixelsDataFrame")) 
        rast <- raster(x)
    else if (missing(coords)) 
        stop("coords must be defined if x is not a SpatialPixelsDataFrame")
    else if (inherits(coords, "SpatialPoints")) {
        if (!isGrid(coords)) 
            stop("this function is only suitable for gridded data")
        coords <- SpatialPoints(coords)
        rast <- raster(SpatialPixelsDataFrame(coords, as.data.frame(x)))
    }
    else {
        if (is.matrix(coords)) 
            coords <- as.data.frame(coords)
        if (is.data.frame(coords)) {
            if (!is.vector(x)) 
                stop("x must be a vector or a SpatialPixelsDataFrame")
            if (!nrow(coords) == length(x)) 
                stop("x and coords must have the same number of elements")
            if (!ncol(coords) == 2) 
                stop("coords must have exactly two columns, for the x and y coordinates of data")
            if (!isGrid(coords)) 
                stop("this function is only suitable for gridded data")
            coords <- SpatialPoints(coords)
            rast <- raster(SpatialPixelDataFrame(coords, as.data.frame(x)))
        }
        else stop("Undefined arguments")
    }
    if (is.null(shape)) 
        plot(rast, zlim = zlim, col = col, ...)
    else {
        plot(shape, col = shapefill, border = shapeborder, ...)
        plot(rast, add = T, zlim = zlim, col = col)
    }
    invisible(rast)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
