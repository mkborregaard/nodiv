\name{nodiv_data}
\alias{nodiv_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nodiv_data(phylo, commatrix, coords, proj4string = CRS(as.character(NA)), type = c("auto", "grid", "points"), shape = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
%%     ~~Describe \code{phylo} here~~
}
  \item{commatrix}{
%%     ~~Describe \code{commatrix} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{proj4string}{
%%     ~~Describe \code{proj4string} here~~
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
function (phylo, commatrix, coords, proj4string = CRS(as.character(NA)), 
    type = c("auto", "grid", "points"), shape = NULL) 
{
    type = match.arg(type)
    if (!(class(phylo) == "phylo")) 
        stop("phylo must be a phylogeny in the ape format")
    if (!class(commatrix) == "distrib_data") {
        if (missing(coords)) 
            stop("if commatrix is not an object of type distrib_data, coords must be specified")
        dist_dat <- distrib_data(commatrix, coords, proj4string, 
            type, shape)
    }
    else dist_dat <- commatrix
    nodiv_dat <- dist_dat
    cat("Comparing taxon names in phylogeny and communities (using picante)\n")
    dat <- match.phylo.comm(phylo, dist_dat$comm)
    nodiv_dat$phylo <- dat$phy
    nodiv_dat$comm <- dat$comm
    if (!(is.data.frame(nodiv_dat$comm) & nrow(nodiv_dat$comm) > 
        1)) 
        stop("The tip labels in the phylogeny do not match the names in the community matrix")
    nodiv_dat$coords <- dist_dat$coords[match(rownames(nodiv_dat$comm), 
        dist_dat$coords$sites), ]
    nodiv_dat$hcom <- matrix2sample(nodiv_dat$comm)
    nodiv_dat$hcom[, 1] <- as.character(nodiv_dat$hcom[, 1])
    nodiv_dat$hcom[, 2] <- as.character(nodiv_dat$hcom[, 2])
    cat("Calculating which species descend from each node\n")
    nodiv_dat$node_species <- Create_node_by_species_matrix(nodiv_dat$phylo)
    class(nodiv_dat) <- c("nodiv_data", "distrib_data")
    return(nodiv_dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
