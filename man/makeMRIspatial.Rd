\name{makeMRIspatial}
\alias{makeMRIspatial}
\title{Obtain Spatial Features of a Mask of an MR Image}
\description{Obtain various spatial features of an MR image, which are
  used in tissue classification.}
\usage{
   makeMRIspatial(mask, nnei, sub, bias)
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
	 the voxel level. The default if \code{FALSE}.}
   \item{bias}{logical; if \code{TRUE}, the spatial parameters for biased
	 field correction are calculated. The default if \code{FALSE}.}
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
  \item{weights}{if \code{bias} is \code{TRUE}, it is a vector
	of weights of neighbors of every voxel for bias field correction;
	otherwise it is equal to \code{NULL}. The default is \code{NULL}.}
  \item{weineighbors}{if \code{bias} is \code{TRUE}, it is a vector
	of sum of weights of neighbors for bias field correction, one
	element per voxel;
	otherwise it is equal to \code{NULL}. The default is \code{NULL}.}
}
\references{
  Dai Feng, Dong Liang, and Luke Tierney (2013)
  An unified Bayesian hierarchical model for MRI tissue classification
  \emph{Statistics in Medicine}

  Dai Feng (2008)
  Bayesian Hidden Markov Normal Mixture Models with Application to MRI
  Tissue Classification
  \emph{Ph. D. Dissertation, The University of Iowa} 
} 
\examples{
  mask <- array(1, dim=c(2,2,2))
  spa <- makeMRIspatial(mask, nnei=6, sub=FALSE)
  spa <- makeMRIspatial(mask, nnei=6, sub=TRUE)
  spa <- makeMRIspatial(mask, nnei=26, sub=TRUE, bias=TRUE)
}
\keyword{spatial}


