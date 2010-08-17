\name{summary.mritc}
\alias{summary.mritc}
\title{Summary Method for Class "mritc"}
\description{
  Summarize some information of an object of class "mritc". 
}
\usage{
   summary.mritc(object, ...)
}
\arguments{
   \item{object}{an object of class "mritc".}
   \item{...}{any additional arguments.}
 }
\value{
  The function computes and returns some
  summary statistics of the object obtained from running the
  function \code{\link[mritc]{mritc}}.
}
\examples{
  T1 <- readMRI(system.file("extdata/t1.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  mask <-readMRI(system.file("extdata/mask.rawb.gz", package="mritc"),
                 c(91,109,91), format="rawb.gz")
  tc.icm <- mritc(T1, mask, method="ICM")
  summary(tc.icm)
}
\keyword{methods}



