\name{parent_pval_representation}
\alias{parent_pval_representation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parent_pval_representation(node_number, rep_matrix, tree = htree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node_number}{
%%     ~~Describe \code{node_number} here~~
}
  \item{rep_matrix}{
%%     ~~Describe \code{rep_matrix} here~~
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
function (node_number, rep_matrix, tree = htree) 
{
    desc = Descendants(node_number)
    if (desc[1] < basal.node(tree) | desc[2] < basal.node(tree)) 
        return(rep(NA, nrow(rep_matrix)))
    desc1row = which(colnames(rep_matrix) == as.character(desc[1]))
    desc2row = which(colnames(rep_matrix) == as.character(desc[2]))
    return(rowMeans(cbind(rep_matrix[, desc1row], 1 - rep_matrix[, 
        desc2row])))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
