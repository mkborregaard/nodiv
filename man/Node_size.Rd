\name{Node_size}
\alias{Node_size}
\alias{Node_sites}
\alias{Node_species}
\alias{Node_occupancy}

\title{
Attributes of individual clades
}
\description{
Functions that summarize the species that descend from a node, and their distribution.
}
\usage{
Node_size(nodiv_data, node = NULL)
Node_sites(nodiv_data, node)
Node_species(nodiv_data, node)
Node_occupancy(nodiv_data, node = NULL)
}

\arguments{
  \item{node}{
A character or numeric variable specifying a certain node in the phylogeny
}
  \item{nodiv_data}{
An object of type \code{nodiv_data} or \code{nodiv_result} 
}
}
\details{
If no \code{node} is specified in \code{Node_size} or \code{Node_occupancy}, the function will calculate the value for all nodes in the phylogeny
}
\value{
For \code{Node_size} and \code{Node_occupancy} the number of species or number of occupied sites of the node; or a vector giving this value for all nodes. For \code{Node_species} or \code{Node_sites} a character vector of names of species descending from a node, or the sites where they occur.
}

\author{
Michael Krabbe Borregaard
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
