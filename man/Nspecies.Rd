\name{Nspecies}
\alias{Nspecies}
\alias{Nsites}
\alias{richness}

\title{
Basic attributes of \code{distrib_data} objects
}
\description{
Calculates the species richness and site occupancy of a \code{distrib_data} object
}
\usage{
Nspecies(x)
Nsites(x)
richness(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of type \code{distrib_data}, \code{nodiv_data} or \code{nodiv_result}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{Nspecies} gives the number of species in the dataset, and \code{Nsites} the number of sites. \code{richness} gives a numeric vector with the species richness of all sites.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Krabbe Borregaard
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line