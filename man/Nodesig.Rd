\name{Nodesig}
\alias{Nodesig}

\title{
Compare the distributions of a pair of sister clades
}
\description{
Compares the overlap of species richness of a pair of sister clades to the expectation from a null model
}
\usage{
Nodesig(nodiv_data, Node_sp, repeats = 100, method = c("quasiswap", "rdtable"), show = F)
}

\arguments{
  \item{nodiv_data}{
An object of type \code{nodiv_data}
}
  \item{Node_sp}{
A boolean vector of the same length as the number of species in nodiv_data, with TRUE for all species belonging to one clade and FALSE for all species belonging to the other
}
  \item{repeats}{
The number of repeats to use for the null model
}
  \item{method}{
The null model algorithm used in the calculations
}
  \item{show}{
Should a progress bar be printed during calculations?
}
}
\details{
The \"rdtable\" algorithm is efficient, but may underestimate the deviation between nodes. It should primarily be used for exploratory data analysis. Note that the more robust \"quasiswap\" algorithm may take hours to run on data sets with a few hundred species.
}
\value{
The return value is a list with 5 elements:
  \item{SR}{A numeric vector of the standardized difference between the empirical and simulated species richness for one node in each site}
  \item{rval}{The rank of the empirical species richness of one node in the distribution of simulated values}
  \item{nodeemp}{The numerical species richness in each site}
  \item{nodemeans}{The mean simulated species richness of one node in each site}
  \item{nodesds}{The standard deviation of the simulated species richness of one node in each site}

}
\references{
Borregaard, M.K., Rahbek, C., Fjeldså, J., Parra, J.L., Whittaker, R.J. & Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12283
}
\author{
Michael Krabbe Borregaard
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Node_analysis}}
}
\examples{
data(Fluvicolines)
desc1 <- Descendants(basal_node(Fluvicolines))[1]
Node_sp <- Fluvicolines$species %in% Node_species(desc1)
res <- Nodesig(Fluvicolines, Node_sp, repeats = 50)
hist(res$SR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
