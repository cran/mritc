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
	dimensional array with the forth dimension equal to 1.}
   \item{file}{the name of the image file to be written out.}
   \item{header}{the header file. For file of the
	 "rawb.gz" format, \code{header} is \code{NULL}.}
   \item{format}{the format of the image file. Right now only the
	 "Analyze", "NIfTI", and raw byte (unsigned with 1 byte per
	 element in the byte stream) gzip formats are supported.}
}
\details{
  Header file is needed for the file of "Analyze" or "NIfTI"
  format. For details of the header files for "Analyze" format,
  see \url{http://eeg.sourceforge.net/ANALYZE75.pdf}.
  For details of the header files for "NIfTI" format and the
  comparison between the "Analyze" and "NIfTI" formates,
  see \url{http://nifti.nimh.nih.gov/}.

  Files of "Analyze" format are written out through
  the function \code{\link[fmri]{write.ANALYZE}}.
  An adapted version of \code{\link[fmri]{write.NIFTI}} is
  used to write files of "NIfTI" format. Instead of ".nii",
  the adapted version writes data out into ".nii.gz" files to save
  space.
}
\seealso{
  \code{\link[fmri]{write.ANALYZE}},
  \code{\link[fmri]{write.NIFTI}}
}
\value{
  Nothing is returned.
}
\examples{
  \dontrun{
  writeMRI(vol, file="vol.rawb.gz", header=NULL, format="rawb.gz")
  
  writeMRI(vol, file="vol", header=niftiheader, format="nifti")

  writeMRI(vol, file="vol", header=analyzeheader, format="analyze")
  }
}
\keyword{utilities}


