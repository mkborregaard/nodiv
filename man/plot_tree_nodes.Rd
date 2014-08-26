\name{plot_tree_nodes}
\alias{plot_tree_nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_tree_nodes(variable, label = variable, tree = htree, main = deparse(substitute(variable)), new.window = FALSE, col = rev(rainbow(256, start = 0.02, end = 0.59)), show.legend = TRUE, sig.cutoff, nodes, roundoff = TRUE, genera.lines = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{
%%     ~~Describe \code{variable} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{new.window}{
%%     ~~Describe \code{new.window} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{show.legend}{
%%     ~~Describe \code{show.legend} here~~
}
  \item{sig.cutoff}{
%%     ~~Describe \code{sig.cutoff} here~~
}
  \item{nodes}{
%%     ~~Describe \code{nodes} here~~
}
  \item{roundoff}{
%%     ~~Describe \code{roundoff} here~~
}
  \item{genera.lines}{
%%     ~~Describe \code{genera.lines} here~~
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
function (variable, label = variable, tree = htree, main = deparse(substitute(variable)), 
    new.window = FALSE, col = rev(rainbow(256, start = 0.02, 
        end = 0.59)), show.legend = TRUE, sig.cutoff, nodes, 
    roundoff = TRUE, genera.lines = FALSE, ...) 
{
    require(ape)
    if (!length(variable) == Nnode(tree)) 
        stop("The length of the variable vector must be the same length as the number of nodes on the tree")
    plotvar <- variable
    if (roundoff & is.numeric(label)) 
        label = round(label, 2)
    if (missing(nodes)) 
        node_index = rep(TRUE, Nnode(tree))
    else {
        node_index <- rep(FALSE, Nnode(tree))
        node_index[nodes] <- TRUE
    }
    node_index[is.na(variable)] <- FALSE
    if (!missing(sig.cutoff)) 
        node_index[variable <= sig.cutoff] <- FALSE
    nodes <- nodenumbers(tree)[node_index]
    plotvar <- plotvar[node_index]
    label <- label[node_index]
    if (new.window) 
        quartz(width = 12, height = 10)
    if (show.legend) {
        oldpar <- par()
        par(mar = c(5, 4, 4, 6) + 0.1)
    }
    plot(tree, show.tip = F, ...)
    title(main)
    nodelabels(text = as.character(label), node = nodes, cex = 0.7, 
        bg = create.cols(plotvar, col))
    if (genera.lines) 
        Put_genera_lines()
    if (show.legend) {
        library(fields)
        par <- oldpar
        image.plot(zlim = range(variable, na.rm = T), col = col, 
            legend.only = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
