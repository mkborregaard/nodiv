\name{toSpatialPoints}
\alias{toSpatialPoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
toSpatialPoints(coords, proj4string, commatrix, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{proj4string}{
%%     ~~Describe \code{proj4string} here~~
}
  \item{commatrix}{
%%     ~~Describe \code{commatrix} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (coords, proj4string, commatrix, type) 
{
    xcol <- 0
    ycol <- 0
    ret <- coords
    colnames(ret) <- tolower(colnames(ret))
    if ("x" \%in\% colnames(ret) & "y" \%in\% colnames(ret)) {
        xcol <- which(colnames(ret) == "x")
        ycol <- which(colnames(ret) == "y")
        ret <- ret[, c(xcol, ycol)]
    }
    else if ("lon" \%in\% substr(colnames(ret), 1, 3) & "lat" \%in\% 
        substr(colnames(ret), 1, 3)) {
        colnames(ret) <- substr(colnames(ret), 1, 3)
        xcol <- which(colnames(ret) == "lon")[1]
        ycol <- which(colnames(ret) == "lat")[1]
        ret <- ret[, c(xcol, ycol)]
    }
    if (!ncol(ret) == 2) 
        stop("ret should be a data.frame or spatial data.frame with 2 columns, giving the x/longitude, and y/latitude of all sites")
    ret <- SpatialPoints(ret, proj4string)
    if (ncol(coords) == 3 & !(xcol + ycol == 0)) 
        sitenames <- coords[, -c(xcol, ycol)]
    else if (nrow(coords) == nrow(commatrix) & !is.null(rownames(commatrix))) 
        sitenames <- rownames(commatrix)
    else if (!is.null(rownames(coords))) 
        sitenames <- rownames(coords)
    else stop("There must be valid site names in the rownames of commatrix or in the coords data")
    type_auto <- ifelse(isGrid(ret), "grid", "points")
    if (type == "auto") 
        type <- type_auto
    else if (!type == type_auto) 
        warning("The specified type of data (points or grid) seems to conflict with the automatic setting. This may cause problems")
    if (type == "grid") 
        ret <- SpatialPixelsDataFrame(ret, data.frame(sites = sitenames, 
            stringsAsFactors = F))
    else ret <- SpatialPointsDataFrame(ret, data.frame(sites = sitenames, 
        stringsAsFactors = F))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
