\name{print.summary_nodiv_data}
\alias{print.summary_nodiv_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.summary_nodiv_data(x, printlen = 4, ...)
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
    cat(paste("Data object with", x$type, "distributions and phylogenetic relationships of", 
        length(x$species), "species in", length(x$sites), "sites\n\n"))
    cat("Species names:\n")
    cat(paste("\t", paste(x$species[1:printlen], collapse = ", "), 
        ", ...\n", sep = ""))
    cat("Site names:\n")
    cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "), 
        ", ...\n\n", sep = ""))
    cat(paste("Species richness:   min:", min(x$richness$richness), 
        "\tmax:", max(x$richness$richness), "\tmean:", round(mean(x$richness$richness), 
            2), "\n"))
    cat(paste("Species occupancy:  min:", min(x$occupancy), "\tmax:", 
        max(x$occupancy), "\tmean:", round(mean(x$occupancy), 
            2), "\n\n"))
    cat(paste("The phylogeny has", length(x$nodes), "internal nodes"))
    if (!is.null(x$node.label)) {
        cat("Node labels:\n")
        if (x$nodes > printlen) {
            cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "), 
                ", ...\n", sep = ""))
        }
        else print(x$node.label)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
