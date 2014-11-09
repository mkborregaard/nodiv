\name{plot_nodes_phylo}
\alias{plot_nodes_phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_nodes_phylo(variable, label = variable, tree, main = deparse(substitute(variable)), zlims, col = brewer.pal(9, "YlOrRd"), show.legend = TRUE, sig.cutoff, nodes, roundoff = TRUE, show.tip.label = F, ...)
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
  \item{zlims}{
%%     ~~Describe \code{zlims} here~~
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
  \item{show.tip.label}{
%%     ~~Describe \code{show.tip.label} here~~
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
function (variable, label = variable, tree, main = deparse(substitute(variable)), 
    zlims, col = brewer.pal(9, "YlOrRd"), show.legend = TRUE, 
    sig.cutoff, nodes, roundoff = TRUE, show.tip.label = F, ...) 
{
    if (!length(variable) == Nnode(tree)) 
        stop("The length of the variable vector must be the same length as the number of nodes on the tree")
    plotvar <- variable
    if (roundoff & is.numeric(label)) 
        label = round(label, 2)
    if (missing(zlims)) 
        zlims <- c(min(plotvar, na.rm = T), max(plotvar, na.rm = T))
    sizes <- par("cex") * 4 * sqrt((plotvar - zlims[1])/zlims[2])
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
    sizes <- sizes[node_index]
    if (show.legend) {
        oldpar <- par()
        par(mar = c(5, 4, 4, 6) + 0.1)
    }
    plot(tree, show.tip.label = show.tip.label, ...)
    title(main)
    nodelabels(pch = 16, node = nodes, col = create.cols(plotvar, 
        col, zlims), cex = sizes)
    nodelabels(text = as.character(label), node = nodes, cex = 0.6, 
        frame = "none")
    if (show.legend) {
        library(fields)
        par <- oldpar
        image.plot(zlim = range(variable, na.rm = T), col = col, 
            legend.only = T, zlim = zlims, smallplot = c(0.85, 
                0.87, 0.38, 0.65))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
