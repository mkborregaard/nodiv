\name{nodiv_res}
\alias{nodiv_res}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nodiv_res(results, nodiv_data, repeats, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
%%     ~~Describe \code{results} here~~
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
function (results, nodiv_data, repeats, method) 
{
    ret <- nodiv_data
    class(ret) <- c("nodiv_result", class(nodiv_data))
    ret$method <- method
    ret$repeats <- repeats
    SR <- sapply(results, "[[", 1)
    ret$SOS <- sapply(nodenumbers(nodiv_data), function(node) parent_representation(node, 
        SR, nodiv_data))
    rval <- sapply(results, "[[", 2)
    par_rval <- sapply(nodenumbers(nodiv_data), function(node) parent_representation(node, 
        rval, nodiv_data))
    pval <- apply(par_rval, 2, pval_sig)
    pval[pval > 1 - 2/repeats] <- 1 - 2/repeats
    pval[pval < 2/repeats] <- 2/repeats
    ret$GND <- apply(pval, 2, function(x) 1 - inv_logit(mean(logit(x), 
        na.rm = T)))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
