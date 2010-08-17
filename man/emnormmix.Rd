\name{emnormmix}
\alias{emnormmix}
\title{Estimate the Parameters of a Normal Mixture Model Using the EM Algorithm }
\description{Use the EM Algorithm to estimate the parameters of a normal mixture model.}
\usage{
   emnormmix(y, prop, mu, sigma, err, maxit, verbose)
}
\arguments{
   \item{y}{vector of observations.}
   \item{prop}{vector of initial estimate of the proportions of
	 different components of a normal mixture model.
     }
	 \item{mu}{vector of initial estimate of the means of different
	   components of a normal mixture model.
     }
   \item{sigma}{vector of initial estimate of the standard deviations of
	 different components of a normal mixture model.}
   \item{err}{relative maximum error(s) used to decide when to stop the
     iteration. It could be a vector of length three corresponding to
     the relative maximum errors of the means, standard deviations, and
     proportions of all components of a normal mixture model. When it is
     a scalar, all have the same relative maximum error.}
   \item{maxit}{maximum number of iterations to perform.}
   \item{verbose}{logical. If \code{TRUE}, then indicate
	 the level of output as the algorithm runs after every 10 iterations.}
}
\details{It is tailor-made for the case when observations are from
  a finite set (MRI data for example). The table lookup method is used to speed up the
  computation.
}
\value{
  \item{prop}{a vector of estimated proportions of different components of a
  normal mixture model.}
 \item{mu}{a vector of estimated means of different components of a
  normal mixture model.}
  \item{sigma}{a vector of estimated standard deviations of different components of a
    normal mixture model.}
}
\examples{
  prop <- c(0.3, 0.3, 0.4)
  mu <- c(-10, 0, 10)
  sigma <- rep(1, 3)
  y<- floor(rnormmix(n=100000, prop, mu, sigma)[,1])
  initial <- initOtsu(y, 2)
  emnormmix(y=y, prop=initial$prop, mu=initial$mu, sigma=initial$sigma,
            err=1e-7, maxit=100, verbose=TRUE)
}
\keyword{utilities}


