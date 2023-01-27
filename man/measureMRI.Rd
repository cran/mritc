\name{measureMRI}
\alias{measureMRI}
\title{Compare the Predicted Classsification Results with the Truth}
\description{Calculate and demonstrate different measures for
  classification results based on the truth.}
\usage{
   measureMRI(intvec, actual, pre)
}
\arguments{
  \item{intvec}{a vector of intensity values. If it is not \code{NULL},
	the density plots of each component corresponding to the actual and
	predicted classification results are shown. The default is
	\code{NULL}.}
   \item{actual}{matrix of the true classification result. Each row
	 corresponds to one voxel. Column \eqn{i} represents the
	 probabilities that all voxels are allocated to tissue type \eqn{i}.}
   \item{pre}{matrix of the predicted classification result. Each row
	 corresponds to one voxel. Column \eqn{i} represents the
	 probabilities that all voxels are allocated to tissue type \eqn{i}.
   }
}
\value{
  \item{mse}{mean square error.}
  \item{misclass}{mis-classification rate.}
  \item{rseVolume}{root square error of volume
	with respect to reference tissue volume.}
  \item{DSM}{Dice Similary Measure of each tissue type.
	\deqn{
  DSM_{a,b}^{t}=\frac{2 \times N_{a \cap b}^t}{N_a^t+N_b^t}  
    }
where \eqn{N_a^t} and \eqn{N_b^t} are the number of voxels classified as
tissue type \eqn{t} by
method \eqn{a} and \eqn{b} respectively, and \eqn{N_{a \cap b}^t} is the number of voxels
classified as tissue type \eqn{t} by both methods \eqn{a} and \eqn{b}.  The larger the
DSM, the more
similar the results from the two methods.
}
  \item{conTable}{confusion table. Each column of the table represents
	the instances in an actual class, while each row represents the
	instances in a predicted class.}
}
\examples{
  #Example 1
  prop <- c(.3, .4, .3)
  mu <- c(-4, 0, 4)
  sigma <- rep(1, 3)
  y <- rnormmix(n=1e4, prop, mu, sigma)
  intvec <- y[,1]
  actual <- y[,2]
  pre <- actual
  pre[sample(1:1e4, 100, replace=FALSE)]  <- sample(1:3, 100, replace=TRUE)
  actual <- do.call(cbind, lapply(1:3, function(i) ifelse(actual==i, 1, 0)))
  pre <- do.call(cbind, lapply(1:3, function(i) ifelse(pre==i, 1, 0)))
  measureMRI(intvec, actual, pre)

   \donttest{
  #Example 2
  T1 <- readMRI(system.file("extdata/t1.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  mask <-readMRI(system.file("extdata/mask.rawb.gz", package="mritc"),
                 c(91,109,91), format="rawb.gz")
  tc.icm <- mritc(T1, mask, method="ICM")

  csf <- readMRI(system.file("extdata/csf.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  gm <- readMRI(system.file("extdata/gm.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  wm <- readMRI(system.file("extdata/wm.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  truth <- cbind(csf[mask==1], gm[mask==1], wm[mask==1])
  truth <- truth/255
  measureMRI(T1[mask==1], truth, tc.icm$prob)
  }
}
\keyword{utilities}


