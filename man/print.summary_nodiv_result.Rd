\name{print.summary_nodiv_result}
\alias{print.summary_nodiv_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.summary_nodiv_result(x, printlen = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{printlen}{
%%     ~~Describe \code{printlen} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, printlen = 4, ...) 
{
    cat(paste("Result of nodiv analysis on", x$type, "data\n"))
    cat(paste("Repeats:", x$repeats, "\n"))
    cat(paste("Null model:", x$method, "\n\n"))
    cat(paste("Species names: (n = ", length(x$species), "):\n", 
        sep = ""))
    cat(paste("\t", paste(x$species[1:printlen], collapse = ", "), 
        ", ...\n\n", sep = ""))
    cat(paste("Site names (n = ", nrow(x$coords), "):\n", sep = ""))
    cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "), 
        ", ...\n\n", sep = ""))
    cat(paste("GND and SOS calculated for", sum(!is.na(x$GND)), 
        "nodes\n"))
    cat(paste("GND values:  min", round(min(x$GND, na.rm = T), 
        2), "max", round(max(x$GND, na.rm = T), 2), "mean", round(mean(x$GND, 
        na.rm = T), 2), "\n\n"))
    cat(paste(length(x$sign), "nodes of", sum(!is.na(x$GND)), 
        "were above 0.6"))
    if (length(x$sign) > 0) {
        cat(":\n")
        for (i in seq_along(x$sign)) {
            nodelab <- ifelse(!is.null(x$phylo$node.label), x$phylo$node.label[x$sign[i]], 
                "")
            cat(paste("\tNode number ", x$nodes[x$sign[i]], ": ", 
                round(x$GND[x$sign[i]], 2), "\t", nodelab, "\n", 
                sep = ""))
        }
    }
    else (cat("\n"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
