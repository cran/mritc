
\name{readMRI}
\alias{readMRI}
\title{Read an MR Image into an Array}
\description{Read an MR image of different formats into an array. 
}
\usage{
   readMRI(file, dim, format)
}
\arguments{
   \item{file}{the name of the image file to be read in.}
   \item{dim}{the dimensions of the image. It is required for the
	 image of type \code{raw.gz}, in which \code{dim} is a vector
	 of length three specifying dimensions in x, y, and z directions.
	 The default is \code{NULL}.}
   \item{format}{the format of the image file. Right now only the
	 "Analyze", "NIfTI", and raw byte (unsigned with 1 byte per
	 element in the byte stream) gzip formats are supported.}
}
\value{
  An array with the appropriate dimensions containing the image
  volume. 
}
\details{
  The files of "Analyze" format are read in through the
  function \code{\link[oro.nifti]{readANALYZE}}.
  The files of "NIfTI" format are read in through
  the function\code{\link[oro.nifti]{readNIfTI}}. 
}
\seealso{
  \code{\link[oro.nifti]{readANALYZE}},
  \code{\link[oro.nifti]{readNIfTI}}
}
\examples{
  \dontrun{
  vol1 <- readMRI("t1.rawb.gz", c(91,109,91), "rawb.gz")
  
  vol2 <- readMRI("t1.nii.gz", format="nifti")

  vol3  <- readMRI("t1.nii", format="nifti")

  vol4 <- readMRI("t1-analyze.img", format="analyze")
  }
}
\keyword{utilities}


