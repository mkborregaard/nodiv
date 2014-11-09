\name{Node_analysis}
\alias{Node_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Node_analysis(nodiv_data, repeats = 100, method = c("rdtable", "quasiswap"), cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodiv_data}{
%%     ~~Describe \code{nodiv_data} here~~
}
  \item{repeats}{
%%     ~~Describe \code{repeats} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
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
function (nodiv_data, repeats = 100, method = c("rdtable", "quasiswap"), 
    cores = 1) 
{
    if (!inherits(nodiv_data, "nodiv_data")) 
        stop("This function only works on objects of class nodiv_data")
    method <- match.arg(method)
    if (cores > 1) {
        stop("The parallel interface is currently unimplemented")
    }
    else {
        pb <- txtProgressBar(min = Nspecies(nodiv_data) + 1, 
            max = Ntip(nodiv_data$phylo) + Nnode(nodiv_data$phylo), 
            style = 3)
        results <- lapply(nodenumbers(nodiv_data), function(node) {
            setTxtProgressBar(pb, node)
            ret <- nodiv_anal(node, nodiv_data, repeats, method)
            ret[, 1:2]
        })
    }
    nodiv_result <- nodiv_res(results, nodiv_data, repeats, method)
    return(nodiv_result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
