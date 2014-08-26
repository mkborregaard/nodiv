\name{plotnode_maps}
\alias{plotnode_maps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotnode_maps(node_number, dispersion, siteresults, new.window = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_number}{
%%     ~~Describe \code{node_number} here~~
}
  \item{dispersion}{
%%     ~~Describe \code{dispersion} here~~
}
  \item{siteresults}{
%%     ~~Describe \code{siteresults} here~~
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
function (node_number, dispersion, siteresults, new.window = FALSE) 
{
    richs <- Node.richness(node_number)
    overrep <- dispersion[, colnames(dispersion) == as.character(node_number)] > 
        0.95
    underrep <- dispersion[, colnames(dispersion) == as.character(node_number)] < 
        0.05
    if (new.window) 
        quartz()
    par(mfrow = c(2, 2))
    map.var(overrep, siteresults$Long, siteresults$Lat, FALSE, 
        main = paste("overrepresented sites for node", node_number), 
        col = c("white", "red2"))
    map.var(underrep, siteresults$Long, siteresults$Lat, FALSE, 
        main = paste("underrepresented sites for node", node_number), 
        col = c("white", "red2"))
    map.var(richs$richness, richs$Long, richs$Lat, FALSE, main = paste("richness of node", 
        node_number))
    map.var(siteresults$richness, siteresults$Long, siteresults$Lat, 
        FALSE, main = "richness of all species")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
