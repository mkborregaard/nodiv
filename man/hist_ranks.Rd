\name{hist_ranks}
\alias{hist_ranks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hist_ranks(node = sample(Nnode(tree), 1) + Ntip(tree), tree = htree, dats = dispersion_corr, comm = hcom, new.window = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{dats}{
%%     ~~Describe \code{dats} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{new.window}{
%%     ~~Describe \code{new.window} here~~
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
function (node = sample(Nnode(tree), 1) + Ntip(tree), tree = htree, 
    dats = dispersion_corr, comm = hcom, new.window = FALSE) 
{
    require(ape)
    nodetree <- extract.clade(tree, node)
    presents <- subset(comm, Species \%in\% nodetree$tip.label)
    presents_index <- which(rownames(dats) \%in\% presents$Plot)
    print(paste("Node", node, "with", length(nodetree$tip.label), 
        "species"))
    if (new.window) 
        quartz()
    hist(dats[, which(as.numeric(colnames(dats)) == node)], main = paste("Node", 
        node, "with", length(nodetree$tip.label), "species"), 
        xlim = c(0, 1), breaks = -2:20/20 + 0.025, xlab = "mean rank in random distribution", 
        col = "green")
    hist(dats[presents_index, which(as.numeric(colnames(dats)) == 
        node)], add = T, breaks = -2:20/20 + 0.025, col = "red")
    abline(v = 0.5, lty = 2)
    abline(v = 0.95, lty = 2)
    abline(v = 0.05, lty = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
