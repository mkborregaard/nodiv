\name{Nspecies}
\alias{Nspecies}
\alias{Nsites}
\alias{richness}
\alias{occupancy}
\alias{sites}
\alias{species}

\title{
Basic attributes of \code{distrib_data} objects
}
\description{
Calculates the species richness and site occupancy of a \code{distrib_data} object
}
\usage{
Nspecies(distrib_data)
Nsites(distrib_data)
richness(distrib_data)
occupancy(distrib_data)
sites(distrib_data)
species(distrib_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distrib_data}{
An object of type \code{distrib_data}, \code{nodiv_data} or \code{nodiv_result}
}
}

\value{
\code{Nspecies} gives the number of species in the dataset, and \code{Nsites} the number of sites. 
\code{richness} gives a numeric vector with the species richness of all sites.
\code{occupancy} gives the number of sites occupied by all species.
\code{sites} gives a list of site names in the object. 
\code{species} gives a list of species names in the object. 
\code{nodes} gives a node labels or node numbers 
}

\author{
Michael Krabbe Borregaard
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
