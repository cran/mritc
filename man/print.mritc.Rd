\name{print.mritc}
\alias{print.mritc}
\title{Print Method for Class "mritc"}
\description{
  Print out some information of an object of class "mritc". 
}
\usage{
   print.mritc(x, ...)
}
\arguments{
   \item{x}{an object of class "mritc".}
   \item{...}{any additional arguments.}
 }
\value{
  The function computes and returns some summary statistics of the
  object obtained from running the function \code{\link[mritc]{mritc}}.
}
\examples{
  T1 <- readMRI(system.file("extdata/t1.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  mask <-readMRI(system.file("extdata/mask.rawb.gz", package="mritc"),
                 c(91,109,91), format="rawb.gz")
  tc.icm <- mritc(T1, mask, method="ICM")
  tc.icm
}
\keyword{methods}



