\name{plottree.dispersion.heights}
\alias{plottree.dispersion.heights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plottree.dispersion.heights(disp = noderesults$overdisp, tree = htree, nodelab = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disp}{
%%     ~~Describe \code{disp} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{nodelab}{
%%     ~~Describe \code{nodelab} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (disp = noderesults$overdisp, tree = htree, nodelab = FALSE, 
    ...) 
{
    nodeages <- max(branching.times(tree)) - branching.times(tree)
    draw.edge <- function(node1, node2) {
        segments(nodeages[node1 - Ntip(tree)], disp[noderesults$nodenumbers == 
            node1], nodeages[node2 - Ntip(tree)], disp[noderesults$nodenumbers == 
            node2])
    }
    plot(nodeages, disp, pch = 16, ...)
    for (i in 1:length(tree$edge[, 1])) {
        if (tree$edge[i, 2] > Ntip(tree)) 
            draw.edge(tree$edge[i, 1], tree$edge[i, 2])
    }
    if (nodelab) 
        text(nodeages - 0.3, disp, 1:Nnode(tree) + Ntip(tree))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
