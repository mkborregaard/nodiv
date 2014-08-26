\name{plot_richness_nodesize}
\alias{plot_richness_nodesize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_richness_nodesize(dats = dispersion, cutoff = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dats}{
%%     ~~Describe \code{dats} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
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
function (dats = dispersion, cutoff = 0.95) 
{
    mean_richness_of_overrepresented_nodes <- numeric()
    for (i in 1:length(siteresults$cell)) {
        site <- siteresults$cell[i]
        represented_nodes <- which(dats[rownames(dats) == as.character(site), 
            ] > cutoff)
        mean_richness_of_overrepresented_nodes[i] <- mean(noderesults$nodesize[represented_nodes])
    }
    plot(siteresults$richness, mean_richness_of_overrepresented_nodes, 
        pch = 16)
    abline(lm(mean_richness_of_overrepresented_nodes ~ siteresults$richness))
    abline(h = mean(noderesults$nodesize), lty = 2)
    print(summary(lm(mean_richness_of_overrepresented_nodes ~ 
        siteresults$richness)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
