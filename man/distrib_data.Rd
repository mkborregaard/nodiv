\name{distrib_data}
\alias{distrib_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
distrib_data(commatrix, coords, proj4string_in = CRS(as.character(NA)), type = c("auto", "grid", "points"), shape = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{commatrix}{
%%     ~~Describe \code{commatrix} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{proj4string_in}{
%%     ~~Describe \code{proj4string_in} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{shape}{
%%     ~~Describe \code{shape} here~~
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
function (commatrix, coords, proj4string_in = CRS(as.character(NA)), 
    type = c("auto", "grid", "points"), shape = NULL) 
{
    type = match.arg(type)
    cat("Checking input data\n")
    if (class(coords) == "SpatialPointsDataFrame" | class(coords) == 
        "SpatialPixelsDataFrame") 
        if (!all.equal(proj4string_in, coords@proj4string)) {
            proj4string_in <- proj4string(coords)
            warning("specified proj4string overridden by the coords data")
        }
    if (is.data.frame(commatrix) & ncol(commatrix) == 3 & !is.numeric(commatrix[, 
        3])) 
        commatrix <- sample2matrix(commatrix)
    if (is.data.frame(commatrix)) 
        commatrix <- as.matrix(commatrix)
    if (!is.matrix(commatrix)) 
        stop("commatrix must be a matrix of 0's and 1's, indicating presence or absence")
    if (!is.numeric(commatrix)) 
        stop("commatrix must be a numeric matrix of 0's and 1's, indicating presence or absence")
    if (!all.equal(sort(unique(as.numeric(commatrix))), 0:1)) 
        stop("commatrix must be a matrix of 0's and 1's, indicating presence or absence")
    if (is.matrix(coords)) 
        coords <- as.data.frame(coords)
    if (is.data.frame(coords)) 
        coords <- toSpatialPoints(coords, proj4string_in, commatrix, 
            type)
    if (class(coords) == "SpatialPixelsDataFrame") 
        type <- "grid"
    else if (class(coords) == "SpatialPointsDataFrame") 
        type <- "points"
    else stop("coords must be a data.frame of coordinates or an sp data.frame object")
    commatrix <- match_commat_coords(commatrix, coords$sites)
    ret <- list(comm = as.data.frame(commatrix), type = type, 
        coords = coords)
    ret$species <- colnames(ret$comm)
    if (!is.null(shape)) 
        ret$shape <- shape
    class(ret) <- "distrib_data"
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
