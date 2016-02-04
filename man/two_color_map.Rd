\name{two_color_map}
\alias{two_color_map}

\title{
Plot the relative amounts of two variables using a purple-green color scheme
}
\description{
Plots site statistics from objects in the \code{nodiv} package
}
\usage{
two_color_map(distrib_data, vec1, vec2, res = 10, showlegend = T, 
legend_coords = c(0.2, 0.26, 0.36, 0.44), type = c("auto", "grid", "points"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distrib_data}{
An object of type \code{distrib_data}, \code{nodiv_data} or \code{nodiv_result}, or, alternatively, an object of spatial points from package \code{sp}, or a matrix of xy values in two columns
}
  \item{vec1}{
A site statistic to be plotted, such as species richness or mean temperature. This can be a vector of values or the name of a site statistic in distrib_data. 
}
  \item{vec2}{
The oter site statistic to be plotted.
}
  \item{res}{
The resolution of colors on the color legend
}
  \item{showlegend}{
A boolean specifying whether to show the legend
}
  \item{legend_coords}{
A vector x1, x2, y1, y2 specifying the location and extent of the color legend
}
  \item{type}{
To specify whether to plot as spatial points or a grid
}
\item{\dots}{
Further arguements to the plotting function
}

}

\value{
A data.frame giving codes and indices of the colors shown in the plot
}

\author{
Michael Krabbe Borregaard
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{distrib_data}}
\code{\link{plot_points}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
