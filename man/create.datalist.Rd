\name{create.datalist}
\alias{create.datalist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create.datalist(dispersion, sitestatistics, coords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dispersion}{
%%     ~~Describe \code{dispersion} here~~
}
  \item{sitestatistics}{
%%     ~~Describe \code{sitestatistics} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
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
function (dispersion, sitestatistics, coords) 
{
    library(ape)
    datalist <- list()
    dispersion[abs(dispersion) == Inf] <- NA
    datalist$siteresults <- merge(sitestatistics, Summarize.sites(dispersion), 
        by = "cell")
    datalist$siteresults <- merge(datalist$siteresults, coords, 
        by = "cell")
    datalist$siteresults <- datalist$siteresults[match(as.character(coords$cell), 
        datalist$siteresults$cell), ]
    datalist$noderesults <- Summarize.nodes(dispersion)
    datalist$rep_matrix <- dispersion[match(as.character(datalist$siteresults$cell), 
        rownames(dispersion)), ]
    datalist$parent_rep_matrix <- parent_representation_matrix(datalist$rep_matrix)
    datalist$parent_rep_matrix[abs(datalist$parent_rep_matrix) == 
        Inf] <- NA
    pod <- apply(datalist$parent_rep_matrix, 2, function(x) mean(abs(x), 
        na.rm = T))
    pod[pod == Inf] <- NA
    datalist$noderesults$parent_overdisp <- pod
    return(datalist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
