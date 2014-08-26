\name{Node.comm}
\alias{Node.comm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Node.comm(node, comm = hcom, tree = htree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
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
function (node, comm = hcom, tree = htree) 
{
    if (node < Ntip(tree)) 
        nodespecs = tree$tip.label[node]
    else nodespecs <- Node.species2(node, tree)
    nodecom <- subset(comm, Species \%in\% nodespecs)
    return(nodecom)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
