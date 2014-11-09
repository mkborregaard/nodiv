\name{Node_occupancy}
\alias{Node_occupancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Node_occupancy(nodiv_data, node = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodiv_data}{
%%     ~~Describe \code{nodiv_data} here~~
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
function (nodiv_data, node = 0) 
{
    if (node == 0) 
        return(sapply(nodenumbers(nodiv_data), function(nod) Node_occupancy(nodiv_data, 
            nod)))
    else if (!node \%in\% nodenumbers(nodiv_data)) 
        stop("node must be one of the internal nodes in nodiv_data$phylo")
    return(length(Node_sites(node, nodiv_data)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
