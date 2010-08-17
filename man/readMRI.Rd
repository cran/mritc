
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
  volume. For the images of "Analyze" format, the ".hdr" file and
  its corresponding ".img" file should be in the same directory.
  For "NIfTI" format, If the extension of the filename is not
  included, the function searches for the "nii.gz" file, then
  for the ".nii" file, and finally for the "hdr/img" files in
  pair.
}
\details{
  The files of "Analyze" format are read in through the
  functions \code{\link[AnalyzeFMRI]{f.read.analyze.volume}}.
  The files of "NIfTI" format are read in through an adapted
  version of \code{\link[fmri]{read.NIFTI}}. Besides the
  original ".nii" files and the "hdr/img" files,
  it can read the "nii.gz" files as well. 
}
\seealso{
  \code{\link[AnalyzeFMRI]{f.read.analyze.volume}},
  \code{\link[fmri]{read.NIFTI}}
}
\examples{
  \dontrun{
  vol1 <- readMRI("t1.rawb.gz", c(91,109,91), "rawb.gz")
  
  vol2 <- readMRI("t1.nii.gz", format="nifti")

  vol3  <- readMRI("t1.nii", format="nifti")

  vol4 <- readMRI("t1-nifti.img", format="nifti")

  vol5 <- readMRI("t1-analyze.img", format="analyze")
  }
}
\keyword{utilities}


