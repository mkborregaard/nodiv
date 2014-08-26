\name{create.cols}
\alias{create.cols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create.cols(vec, col = rev(rainbow(256, start = 0.02, end = 0.59)), zlims = c(min(vec, na.rm = T), max(vec, na.rm = T)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
%%     ~~Describe \code{vec} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{zlims}{
%%     ~~Describe \code{zlims} here~~
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
function (vec, col = rev(rainbow(256, start = 0.02, end = 0.59)), 
    zlims = c(min(vec, na.rm = T), max(vec, na.rm = T))) 
{
    vec = vec - zlims[1]
    vec = floor(vec * (length(col) - 1)/(zlims[2] - zlims[1])) + 
        1
    return(col[vec])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
