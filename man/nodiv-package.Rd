\name{nodiv-package}
\alias{nodiv-package}
\alias{nodiv}
\docType{package}
\title{
nodiv - Node-based analysis of species distributions
}
\description{
The package implements Borregaard et al. (2014) method for identifying nodes in a phylogeny associated with divergent distributions. The main algorithm goes through each node in the phylogeny and relates node overlap to a null model. The package also provides functions for preparing the data sets, for exploratory plots and further analysis
}
\details{
\tabular{ll}{
Package: \tab nodiv\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2014-11-09\cr
License: \tab GPL-2\cr
}
~~ An overview of how to use the package, including the most ~~
~~ important functions ~~
}
\author{
Michael Krabbe Borregaard

Maintainer: Michael Krabbe Borregaard <mkborregaard@snm.ku.dk>
}
\references{
Borregaard, M.K., Rahbek, C., Fjeldså, J., Parra, J.L., Whittaker, R.J. & Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12283
}
~~ Optionally other standard keywords, one per line, from file ~~
~~ KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
}
\examples{
data(Fluvicolines)
res <- Node_analysis(Fluvicolines, 50, "rdtable")
plot(res)
plotSOS(res, 398)
}
