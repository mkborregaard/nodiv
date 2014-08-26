\name{map.var}
\alias{map.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
map.var(variable, Long, Lat, new.window = FALSE, xrange = range(Long), yrange = range(Lat), bordercolor = "darkgrey", cut.to.coast = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{
%%     ~~Describe \code{variable} here~~
}
  \item{Long}{
%%     ~~Describe \code{Long} here~~
}
  \item{Lat}{
%%     ~~Describe \code{Lat} here~~
}
  \item{new.window}{
%%     ~~Describe \code{new.window} here~~
}
  \item{xrange}{
%%     ~~Describe \code{xrange} here~~
}
  \item{yrange}{
%%     ~~Describe \code{yrange} here~~
}
  \item{bordercolor}{
%%     ~~Describe \code{bordercolor} here~~
}
  \item{cut.to.coast}{
%%     ~~Describe \code{cut.to.coast} here~~
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
function (variable, Long, Lat, new.window = FALSE, xrange = range(Long), 
    yrange = range(Lat), bordercolor = "darkgrey", cut.to.coast = FALSE, 
    ...) 
{
    require(fields)
    require(maptools)
    variable[is.na(variable)] <- -9999
    cellsize_x <- min(diff(sort(unique(Long))), na.rm = T)
    cellsize_y <- min(diff(sort(unique(Lat))), na.rm = T)
    map <- as.image(variable, x = cbind(Long, Lat), grid = list(x = seq(xrange[1], 
        xrange[2], by = cellsize_x), y = seq(yrange[1], yrange[2], 
        by = cellsize_y)))
    map$z[map$z == -9999] <- NA
    if (new.window) 
        quartz()
    image.plot(map, ...)
    sa <- readShapePoly("~/shapecoast/World_coast.shp")
    plot(sa, add = T, border = bordercolor)
    if (cut.to.coast) {
        sc <- readShapePoly("~/shapecoast/World_coast.shp")
        plot(sc, add = T, border = NULL, col = "white", usePolypath = TRUE)
        box()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
