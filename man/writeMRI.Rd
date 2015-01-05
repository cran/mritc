\name{writeMRI}
\alias{writeMRI}
\title{Write an MR Image}
\description{Write an MR image into a file of different formats. 
}
\usage{
   writeMRI(data, file, header, format)
}
\arguments{
  \item{data}{MRI data in a three dimensional array or four
	dimensional array with the forth dimension equal to 1.
	It could be also an object of class "nifti" or "anlz" as
	defined in the package
	\code{oro.nifti}.}
   \item{file}{the name of the image file to be written out.}
   \item{header}{the header file. \code{header} is set as \code{NULL}
	 for backward compatibility.}
   \item{format}{the format of the image file. Right now only the
	 "Analyze", "NIfTI", and raw byte (unsigned with 1 byte per
	 element in the byte stream) gzip formats are supported.}
}
\details{
  Header file is not needed for the file of "Analyze" or "NIfTI"
  format anymore. 

  Files of "Analyze" format are written out through
  the function \code{\link[oro.nifti]{writeANALYZE}}.
  Files of "NIfTI" format are written out through
  the function \code{\link[oro.nifti]{writeNIfTI}}.
}
\seealso{
  \code{\link[oro.nifti]{writeANALYZE}},
  \code{\link[oro.nifti]{writeNIfTI}}
}
\value{
  Nothing is returned.
}
\examples{
  \dontrun{
  writeMRI(vol, file="vol.rawb.gz", header=NULL, format="rawb.gz")
  
  writeMRI(vol, file="vol", header=NULL, format="nifti")

  writeMRI(vol, file="vol", header=NULL, format="analyze")
  }
}
\keyword{utilities}


