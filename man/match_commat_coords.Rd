\name{match_commat_coords}
\alias{match_commat_coords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
match_commat_coords(commatrix, sitenames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{commatrix}{
%%     ~~Describe \code{commatrix} here~~
}
  \item{sitenames}{
%%     ~~Describe \code{sitenames} here~~
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
function (commatrix, sitenames) 
{
    cat("Comparing sites in community data and spatial points\n")
    if (is.null(rownames(commatrix))) 
        if (nrow(commatrix) == length(sitenames)) 
            rownames(commatrix) <- sitenames
        else stop("The number of sites in coords and the data matrix do not fit and there are no rownames in the community matrix to use for matching")
    if (sum(sitenames \%in\% rownames(commatrix)) < 2) 
        stop("the coordinate names and the rownammes of the community matrix do not match")
    sitenames <- sitenames[sitenames \%in\% rownames(commatrix)]
    commatrix <- commatrix[match(sitenames, rownames(commatrix)), 
        ]
    return(commatrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
