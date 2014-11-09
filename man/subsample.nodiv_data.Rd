\name{subsample.nodiv_data}
\alias{subsample.nodiv_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subsample.nodiv_data(x, sites = NULL, species = NULL, node = NULL)
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
  \item{node}{
%%     ~~Describe \code{node} here~~
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
function (x, sites = NULL, species = NULL, node = NULL) 
{
    ret_phylo <- x$phylo
    ret_phylo$node.label <- nodenumbers(x)
    if (!node == 0) {
        ret_phylo <- extract.clade(ret_phylo, node)
        species <- ret_phylo$tip.label
    }
    ret <- subsample.distrib_data(x, sites, species)
    new_phylo <- match.phylo.comm(ret_phylo, ret$comm)$phy
    ret$phylo <- drop.tip(x$phylo, which(!x$species \%in\% new_phylo$tip.label))
    ret$hcom <- subset(x$hcom, plot \%in\% ret$coords$sites & id \%in\% 
        ret$species)
    ret$node_species <- x$node_species[, colnames(x$node_species) \%in\% 
        ret$species]
    ret$node_species <- ret$node_species[rowSums(ret$node_species) > 
        0, ]
    ret$node_species <- ret$node_species[as.numeric(rownames(ret$node_species)) >= 
        as.numeric(new_phylo$node.label)[1], ]
    rownames(ret$node_species) <- nodenumbers(ret$phylo)
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
