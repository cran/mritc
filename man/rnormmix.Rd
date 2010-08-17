\name{rnormmix}
\alias{rnormmix}
\title{Generate Random Samples from a Normal Mixture Model}
\description{Generate random samples from a normal mixture model.}
\usage{
   rnormmix(n, prop, mu, sigma)
}
\arguments{
   \item{n}{number of observations.}
   \item{prop}{a vector of proportions of different components.}
   \item{mu}{a vector of means of different components.}
   \item{sigma}{a vector of standard deviations of different components}
}
\value{
  A matrix with each row corresponding to one sample. The first column are
  sample values from a normal mixture model; the second column are the
  components from which observations come.
}
\examples{
  prop <- c(.17, .48, .35)
  mu <- c(-4, 0, 4)
  sigma <- rep(1, 3)
  y <- rnormmix(n=10000, prop, mu, sigma)
  densityplot(~ y[,1], groups = y[,2],
              plot.points = FALSE, ref = TRUE,
              xlab="sample values",
              auto.key = list(columns = 3))
}
\keyword{distribution}


