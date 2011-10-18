\name{plot.mritc}
\alias{plot.mritc}
\title{Plot Method for Class "mritc"}
\description{
  Visualize MRI tissue classification results. 
}
\usage{
   \method{plot}{mritc}(x, ...)
}
\arguments{
   \item{x}{an object of class "mritc"}
   \item{...}{any additional arguments for function \code{\link[misc3d]{slices3d}}.}
}
\details{
  Allocate a voxel to the tissue type with the highest probability and
  then use \code{\link[misc3d]{slices3d}} to show the result.
}
 \value{
   \code{NULL}
}
\seealso{
  \code{\link[misc3d]{slices3d}}
}
\examples{
  T1 <- readMRI(system.file("extdata/t1.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  mask <-readMRI(system.file("extdata/mask.rawb.gz", package="mritc"),
                 c(91,109,91), format="rawb.gz")
  tc.icm <- mritc(T1, mask, method="ICM")
  plot(tc.icm)
}
\keyword{methods}



