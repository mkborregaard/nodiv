\name{plotnode_maps_new_parent}
\alias{plotnode_maps_new_parent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotnode_maps_new_parent(node_number, par_rep_matrix = parent_rep_matrix, tree = htree, comm = hcom, coords = dat.LL, new.window = FALSE, plottype = c("map", "points"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_number}{
%%     ~~Describe \code{node_number} here~~
}
  \item{par_rep_matrix}{
%%     ~~Describe \code{par_rep_matrix} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{new.window}{
%%     ~~Describe \code{new.window} here~~
}
  \item{plottype}{
%%     ~~Describe \code{plottype} here~~
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
function (node_number, par_rep_matrix = parent_rep_matrix, tree = htree, 
    comm = hcom, coords = dat.LL, new.window = FALSE, plottype = c("map", 
        "points")) 
{
    require(ape)
    if (new.window) 
        quartz()
    par(mfrow = c(2, 2))
    plotoverrep(node_number, par_rep_matrix, coords, tree, new.window = FALSE, 
        plottype = plottype)
    plotnode(node_number, comm, tree, coords, new.window = FALSE, 
        plottype = plottype, main = paste("richness of node", 
            node_number))
    plotnode(Descendants(node_number)[1], comm, tree, coords, 
        new.window = FALSE, plottype = plottype, main = paste("richness of descendant 1", 
            Descendants(node_number)[1]))
    plotnode(Descendants(node_number)[2], comm, tree, coords, 
        new.window = FALSE, plottype = plottype, main = paste("richness of descendant 2", 
            Descendants(node_number)[2]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
