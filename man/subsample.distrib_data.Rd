\name{subsample.distrib_data}
\alias{subsample.distrib_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subsample.distrib_data(x, sites = NULL, species = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sites}{
%%     ~~Describe \code{sites} here~~
}
  \item{species}{
%%     ~~Describe \code{species} here~~
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
function (x, sites = NULL, species = NULL) 
{
    if (inherits(sites, "SpatialPoints")) 
        sites <- as.character(sites@data)
    keep_sites <- F
    if (is.character(sites)) {
        if (length(sites) == 1) {
            if (sites == "all") 
                keep_sites <- T
        }
        else sites <- match(sites, x$coords$sites)
    }
    if (is.null(sites) | keep_sites) 
        sites <- 1:Nsites(x)
    keep_species <- F
    if (length(species) == 1) {
        if (species == "all") 
            keep_species <- T
    }
    else species <- match(species, x$species)
    if (is.null(species) | keep_species) 
        species <- 1:Nspecies(x)
    ret <- x
    ret$comm <- ret$comm[sites, species]
    if (keep_sites) 
        sites_keep <- rep_len(TRUE, nrow(ret$comm))
    else sites_keep <- (rowSums(ret$comm, na.rm = T) > 0)
    if (keep_species) 
        sites_keep <- rep_len(TRUE, nrow(ret$comm))
    else species_keep <- (colSums(ret$comm, na.rm = T) > 0)
    ret$comm <- ret$comm[sites_keep, species_keep]
    ret$species <- colnames(ret$comm)
    ret$coords <- ret$coords[ret$coords$sites \%in\% rownames(ret$comm), 
        ]
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
