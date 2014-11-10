\name{basal_node}
\alias{basal_node}
\alias{nodenumbers}
\alias{Parent}
\alias{Descendants}
\alias{Sister}

\title{
Identify nodes in a phylogeny
}
\description{
Returns the internal (\code{ape}) node numbers corresponding to certain nodes in the phylogeny
}
\usage{
basal_node(tree)
nodenumbers(tree)
Descendants(node, tree)
Parent(node, tree)
Sister(node, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
An object of type \code{phylo}, or alternatively an object of type \code{nodiv_data}
}
  \item{node}{
The internal node number of a focal node in the phylogeny
}
}
\details{
\code{Parent} gives the code immediately basal to the focal node, \code{Descendants} gives the nodes descending from the focal node, and \code{Sister} gives the sister node(s).
}
\value{
The internal (\code{ape}) node number of the node
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
