\name{Sitestats}
\alias{Sitestats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sitestats(comm = hcom, tree = htree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
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
function (comm = hcom, tree = htree) 
{
    require(ape)
    cell <- sort(unique(comm$Plot))
    richness = numeric()
    mean_range_in_cell = numeric()
    ranges = table(comm$Species)
    ranges = ranges[match(tree$tip.label, names(ranges))]
    for (plotnum in seq_along(cell)) {
        print(plotnum)
        testcom <- comm[comm$Plot == cell[plotnum], ]
        emp <- as.character(testcom$Species)
        mean_range_in_cell[plotnum] = mean(ranges[tree$tip.label \%in\% 
            emp])
        richness[plotnum] = length(emp)
    }
    return(data.frame(cell, richness, mean_range_in_cell))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
