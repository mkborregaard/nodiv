\name{plotsite}
\alias{plotsite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotsite(site, dats = dispersion_corr, tree = htree, comm = hcom, genera_lines = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{site}{
%%     ~~Describe \code{site} here~~
}
  \item{dats}{
%%     ~~Describe \code{dats} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{genera_lines}{
%%     ~~Describe \code{genera_lines} here~~
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
function (site, dats = dispersion_corr, tree = htree, comm = hcom, 
    genera_lines = TRUE) 
{
    require(ape)
    quartz(width = 20, height = 12)
    layout(matrix(c(1, 2), ncol = 2), width = c(4, 4))
    spec_index = which(tree$tip.label \%in\% Site.species(site, 
        comm))
    colortip <- rep("transparent", Ntip(tree))
    colortip[spec_index] <- "black"
    plot(tree, show.tip = F)
    tiplabels(tip = spec_index, frame = "none", pch = 16, col = "black", 
        bg = colortip, adj = 1, cex = 1)
    overreps = nodenumbers(tree)[dats[site, ] > 1.96]
    underreps = nodenumbers(tree)[dats[site, ] < -1.96]
    if (length(overreps) > 0) 
        nodelabels(node = overreps, cex = 0.9, text = " + ", 
            font = 2, bg = "black", frame = "rect", adj = c(1.1, 
                0), col = "white")
    if (length(underreps) > 0) 
        nodelabels(node = underreps, cex = 0.9, text = " - ", 
            font = 2, bg = "white", frame = "rect", adj = c(-0.2, 
                0), col = "black")
    if (genera_lines) 
        Put_genera_lines()
    plot(siteresults$Long, siteresults$Lat, type = "n")
    background()
    points(siteresults$Long[site], siteresults$Lat[site], col = "red", 
        pch = 16, cex = 1.5)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
