\name{nodiv_anal}
\alias{nodiv_anal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nodiv_anal(node, nodiv_data, repeats, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{nodiv_data}{
%%     ~~Describe \code{nodiv_data} here~~
}
  \item{repeats}{
%%     ~~Describe \code{repeats} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (node, nodiv_data, repeats, method) 
{
    if (node == basal_node(nodiv_data)) 
        return(data.frame(ses = rep(NA, Nsites(nodiv_data)), 
            pval = rep(NA, Nsites(nodiv_data)), nodeemp = rep(NA, 
                Nsites(nodiv_data)), nodemeans = rep(NA, Nsites(nodiv_data)), 
            nodesds = rep(NA, Nsites(nodiv_data))))
    parent_data <- subsample(nodiv_data, node = Parent(node, 
        nodiv_data))
    parNode_sites <- match(parent_data$coords$sites, nodiv_data$coords$sites)
    Node_sp <- (parent_data$species \%in\% Node_species(node, nodiv_data))
    res_object <- Nodesig(parent_data, Node_sp, repeats, method)
    res_object <- lapply(res_object, function(x) {
        ret <- rep(NA, Nsites(nodiv_data))
        ret[parNode_sites] <- x
        ret
    })
    return(as.data.frame(res_object))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
