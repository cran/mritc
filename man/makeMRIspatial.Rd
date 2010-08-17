\name{makeMRIspatial}
\alias{makeMRIspatial}
\title{Obtain Spatial Features of a Mask of an MR Image}
\description{Obtain various spatial features of an MR image, which are
  used in tissue classification.}
\usage{
   makeMRIspatial(mask, nnei, sub)
}
\arguments{
   \item{mask}{three dimensional array. The voxels with value 1 are
     inside the mask; with value 0 are outside. We just focus on voxels
inside the mask.}
   \item{nnei}{the number of neighbors. Right now only 6, 18,
	 and 26 neighbors are supported. For a 3D image,
	 besides defining 6 neighbors in the x, y, and z directions,
	 one can add 12 diagonal neighbors in the x-y, x-z, and y-z
	 planes, and another 8 on the 3D diagonals. This leads to a
	 six neighbor structure, an eighteen neighbor structure, and
	 a twenty-six neighbor structure.}
   \item{sub}{logical; if \code{TRUE}, a new \code{mask}
	 which splits each voxel into
	 eight subvoxels is generated, and then obtain the neighbors and blocks
	 of subvoxels; otherwise obtain the neighbors and blocks at
	 the voxel level.}
}
\value{
  A list containing the following components:
  \item{neighbors}{a matrix, each row of which giving the neighbors of a
    voxel or subvoxel. The number of rows is equal to the number of
	(sub)voxels within the \code{mask} and the number of columns is the
	number of neighbors of each (sub)voxel. For the (sub)voxels on the
	boundaries, when one or more of their
	neighbors are missing, the missing are represented by the total
	number of (sub)voxels within the \code{mask} plus 1.}	
  \item{blocks}{the (sub)voxels within each block are mutually
	independent given the (sub)voxels in other blocks.}
  \item{sub}{logical; the same as the input \code{sub}.}
  \item{subvox}{if \code{sub} is \code{TRUE}, it is a matrix,
	with each row giving the eight subvoxels of a voxel;
	otherwise it is equal to \code{NULL}.}
}
\references{
 Dai Feng (2008)
 Bayesian Hidden Markov Normal Mixture Models with Application to MRI
 Tissue Classification
 \emph{Ph. D. Dissertation, The University of Iowa} 
} 
\examples{
  mask <- array(1, dim=c(2,2,2))
  spa <- makeMRIspatial(mask, nnei=6, sub=FALSE)
  spa <- makeMRIspatial(mask, nnei=6, sub=TRUE)
}
\keyword{spatial}


