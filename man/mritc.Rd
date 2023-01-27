\name{mritc}
\title{MRI Tissue Classification Using Various Methods}
\alias{mritc.em}
\alias{mritc.icm}
\alias{mritc.hmrfem}
\alias{mritc.pvhmrfem}
\alias{mritc.bayes}
\alias{mritc}
\description{
  Conduct the MRI tissue classification using different methods including:
  the normal mixture model (NMM) fitted by the Expectation-Maximization (EM)
  algorithm;
  the hidden Markov normal mixture model (HMNMM) fitted by the Iterated
  Conditional Mode (ICM) algorithm, the Hidden Markov
  Random Field EM (HMRFEM) algorithm, or the Bayesian Markov chain Monte
  Carlo method (MCMC); 
  the partial volume HMNMM fitted by the modified EM (PVHMRFEM)
  algorithm or the higher resolution HMNMM fitted by the MCMC method (MCMCsub);
  the HMNMM with both PV and intensity non-uniformity addressed (MCMCsubbias).
}
\usage{
   mritc.em(y, prop, mu, sigma, err, maxit, verbose)
   mritc.icm(y, neighbors, blocks, spatialMat, beta, mu, sigma,
             err, maxit, verbose)
   mritc.hmrfem(y, neighbors, blocks, spatialMat, beta, mu, sigma,
                err, maxit, verbose)
   mritc.pvhmrfem(y, neighbors, blocks, spatialMat, beta, mu, sigma,
              err, maxit, verbose)
   mritc.bayes(y, neighbors, blocks, sub, subvox,
               subbias, neighbors.bias, blocks.bias, weineighbors.bias, weights.bias,
               spatialMat, beta, mu, sigma, niter, verbose)
   mritc(intarr, mask, method)
}
\arguments{
  \item{y}{a vector of intensity values of voxels.}
  \item{prop}{a vector of initial estimate of the proportions of
	different components of a normal mixture model. It can be obtained using
    the function \code{\link[mritc]{initOtsu}}.}
  \item{mu}{a vector of initial estimate of the means of different
	components of a normal mixture model. It can be obtained using
    the function \code{\link[mritc]{initOtsu}}.}
  \item{sigma}{a vector of initial estimates of the standard
	deviations of different components of a normal mixture model.
	It can be obtained using the function \code{\link[mritc]{initOtsu}}.}
  \item{err}{relative maximum error(s) used to decide when to stop the
    iteration. It could be a vector corresponding to
    the relative maximum errors of the means, standard deviations (for
	\code{\link[mritc]{mritc.em}}, \code{\link[mritc]{mritc.icm}},
	\code{\link[mritc]{mritc.hmrfem}},
	and \code{\link[mritc]{mritc.pvhmrfem}}), and
    proportions (for \code{\link[mritc]{mritc.em}})
	of all components of a normal mixture model. When it is
    a scalar, all have the same relative maximum error. The default
	value is 1e-4.
  }
  \item{maxit}{maximum number of iterations to perform. The default is
	200 for \code{\link[mritc]{mritc.em}}, 20 for
	\code{\link[mritc]{mritc.icm}}, \code{\link[mritc]{mritc.hmrfem}},
	and \code{\link[mritc]{mritc.pvhmrfem}}.} 
  \item{verbose}{logical. If \code{TRUE}, then indicate
	the level of output as the algorithm runs.}
  \item{neighbors}{a matrix of neighbors of voxels. One row per voxel.
	It can be obtained 
	using the function \code{\link[mritc]{makeMRIspatial}}.}
  \item{blocks}{split voxels into different blocks to use the
    checker-board idea. It can be obtained using the function
	\code{\link[mritc]{makeMRIspatial}}.}
  \item{spatialMat}{a matrix defining the spatial relationship
    in a Potts model. The default value is \code{diag(1,3)}
	for three components models for \code{\link[mritc]{mritc.icm}},
	\code{\link[mritc]{mritc.hmrfem}} and
	\code{\link[mritc]{mritc.bayes}} when \code{sub} is \code{FALSE}
	and \code{matrix(c(2,0,-1,0,2,0,-1,0,2), nrow=3)}
	when \code{sub} is \code{TRUE}.
	For \code{\link[mritc]{mritc.pvhmrfem}} the default is
	\code{matrix(c(2, 1, -1, -1, -1, 1, 2, 1, -1, -1, -1, 1,
	  2, 1, -1, -1, -1, 1, 2, 1, -1, -1, -1, 1, 2), ncol=5)}.}
  \item{beta}{the parameter 'inverse temperature' of the Potts
	model. The default value is 0.4 for \code{\link[mritc]{mritc.icm}},
	0.5 for \code{\link[mritc]{mritc.hmrfem}}, 0.6 for
	\code{\link[mritc]{mritc.pvhmrfem}}.
	For \code{\link[mritc]{mritc.bayes}}, the default is 0.7
	when \code{sub} is \code{FALSE} and 0.3 when \code{sub} is
	\code{TRUE}.}
  \item{sub}{logical; if \code{TRUE}, use the higher resolution model;
	otherwise, use the whole voxel method.}
  \item{subvox}{for \code{\link[mritc]{mritc.bayes}},
	the match up tabel of voxels and their corresponding
     subvoxels for the higher resolution model. It can be
	 obtained using the function \code{\link[mritc]{makeMRIspatial}}.
	 For the whole voxel method, \code{subvox}=\code{NULL}}
  \item{subbias}{logical; if \code{TRUE}, use the model that addresses
	both the PV and intensity non-uniformity. The default is \code{FALSE}.}
  \item{neighbors.bias}{a matrix of neighbors of bias field. One row per voxel.
	It can be obtained using the function
	\code{\link[mritc]{makeMRIspatial}}. The default is \code{NULL}.}
  \item{blocks.bias}{blocks for bias field. It can be obtained using the function
	\code{\link[mritc]{makeMRIspatial}}. The default is \code{NULL}.}
  \item{weineighbors.bias}{a vector of sum of weights of neighbors of bias
	field. One element per voxel. It can be obtained using the function
	\code{\link[mritc]{makeMRIspatial}}. The default is \code{NULL}.}
  \item{weights.bias}{a vector of weights of different neighbors of
	every voxel. It can be obtained using the function
	\code{\link[mritc]{makeMRIspatial}}. The default is \code{NULL}.}
  \item{niter}{the number of iterations for
	 \code{\link[mritc]{mritc.bayes}}. The default values are 1000 and
	 100 for with and without bias field correction, respectively.
	 The default values seem to be adequate in many cases.}
  \item{intarr}{a three dimensional array of an MR image.}
  \item{mask}{a mask of the MR image. Voxels with value 1 are
	inside the	brain and value 0 are outside. Focus on voxels
	within the brain.}
  \item{method}{a string giving the method for MRI tissue
	classification. It must be one of "EM", "ICM", "HMRFEM", 
	"MCMC", "PVHMRFEM", "MCMCsub", or "MCMCsubbias" corresponding to using the
	NMM fitted by the EM algorithm;
    the HMNMM fitted by the ICM algorithm, the HMRFEM algorithm,
	or the MCMC; the partial volume HMNMM fitted by the PVHMRFEM algorithm;
    the higher resolution HMNMM fitted by the MCMC; the HMNMN addressing
	both the PV and intensity non-uniformity.
    It can be abbreviated. The default is "EM".
  }
}
\details{
  The function \code{\link[mritc]{mritc}} integrates functions
  \code{\link[mritc]{mritc.em}},
  \code{\link[mritc]{mritc.icm}}, \code{\link[mritc]{mritc.hmrfem}}, 
  \code{\link[mritc]{mritc.pvhmrfem}}, and \code{\link[mritc]{mritc.bayes}}.
  It provides a uniform platform with easier
  usage. The user just need to specify the input MR image, the
  \code{mask} of the image, and the \code{method} used.
  The other parameters are specified automatically as follows.
  The parameters for the initial estimates of the proportions,
  means, and standard deviations of the normal mixture model are
  obtained using the function \code{\link[mritc]{initOtsu}}.
  As to the parameters related to the Potts model,
  the six neighbor structure is used and then the \code{neighbors},
  \code{blocks}, and \code{subvox} are
  obtained using the function \code{\link[mritc]{makeMRIspatial}}.
  For the bias field correction, the twenty-six neighbor structure is
  used and then the \code{neighbors.bias}, \code{blocks.bias},
  \code{weineighbors.bias} and \code{weights.bias} are
  obtained using the function \code{\link[mritc]{makeMRIspatial}}.
  The other parameters are taken as the default values
  for each method. The process is reported during iterations.
  
}
\note{
  The functions support the various normal mixture models with
  at most eight components.
 }
\value{
  For \code{\link[mritc]{mritc}}, it generates
  an object of class "mritc" which is a list containing the
  following components:
   \item{prob}{a matrix, one row per voxel and each column
	corresponding to the probabilities of being allocated to each
	component of a normal mixture model.}
  \item{mu}{a vector of estimated means of the normal mixture model.}
  \item{sigma}{a vector of estimated standard deviations of the
	 normal mixture model.}
   \item{method}{the method used for computation.}
   \item{mask}{mask of an brain. Voxels inside it are classified.}
   
   Generic functions \code{\link[mritc]{print.mritc}}, 
   \code{\link[mritc]{summary.mritc}}, and
   \code{\link[mritc]{plot.mritc}} are provided.

   For others, only \code{prob}, \code{mu}, and \code{sigma} are
   generated.
 }
\references{
  Julian Besag (1986)
  On the statistical analysis of dirty pictures (with discussion)
  \emph{Journal of the Royal Statistical Society. Series B
	(Methodological)}
  \bold{vol. 48} 259-302
  
  Meritxell Bach Cuadra, Leila Cammoun, Torsten Butz, Olivier Cuisenaire,
  and Jean-Philippe Thiran (2005)
  Comparison and validation of tissue modelization and statistical
  classification methods in T1-weighted {MR} brain images
  \emph{IEEE Transactions on Medical Imaging}, \bold{vol.24} 1548-1565

  Dai Feng, Dong Liang, and Luke Tierney (2014)
  An unified Bayesian hierarchical model for MRI tissue classification
  \emph{Statistics in Medicine}
  \bold{vol.33, issue 8} 1349-1368
  
  Dai Feng, Luke Tierney, and Vincent Magnotta (2012)
  MRI tissue classification using high resolution
  Bayesian hidden Markov normal mixture models
  \emph{Journal of the American Statistical Association},
  \bold{vol.107, no.497} 102-119

  Dai Feng and Luke Tierney (2011)
  mritc: A package for MRI tissue classification
  \emph{Journal of Statistical Software}, \bold{vol.44, no.7} 1-20
  \url{https://www.jstatsoft.org/v44/i07/}

  Dai Feng (2008)
  Bayesian hidden Markov normal mixture models with application to MRI
  tissue classification
  \emph{Ph. D. Dissertation, The University of Iowa}

  Yongyue Zhang, Michael Brady, and Stephen Smith (2001)
  Segmentation of brain MR images through a hidden Markov random field
  model and the expectation-maximization algorithm
  \emph{IEEE Transactions on Medical Imaging} \bold{vol. 20} 45-57
 
} 
\examples{
  #Example 1
  T1 <- readMRI(system.file("extdata/t1.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  mask <- readMRI(system.file("extdata/mask.rawb.gz", package="mritc"),
                  c(91,109,91), format="rawb.gz")
  y <- T1[mask==1]
  initial <- initOtsu(y, 2)
  prop <- initial$prop
  mu <- initial$mu
  sigma <- initial$sigma
  tc.em <- mritc.em(y, prop, mu, sigma, verbose=TRUE)
 
  \donttest{
  mrispatial <- makeMRIspatial(mask, nnei=6, sub=FALSE)
  tc.icm <- mritc.icm(y, mrispatial$neighbors, mrispatial$blocks,
                      mu=mu, sigma=sigma, verbose=TRUE)
  tc.hmrfem <- mritc.hmrfem(y, mrispatial$neighbors, mrispatial$blocks,
                            mu=mu, sigma=sigma, verbose=TRUE)
  tc.pvhmrfem <- mritc.pvhmrfem(y, mrispatial$neighbors, mrispatial$blocks,
                                mu=mu, sigma=sigma, verbose=TRUE)
  tc.mcmc <- mritc.bayes(y, mrispatial$neighbors, mrispatial$blocks,
                         mrispatial$sub, mrispatial$subvox,
                         mu=mu, sigma=sigma, verbose=TRUE)

  mrispatial <- makeMRIspatial(mask, nnei=6, sub=TRUE)
  tc.mcmcsub <- mritc.bayes(y, mrispatial$neighbors, mrispatial$blocks,
                         mrispatial$sub, mrispatial$subvox,
                         mu=mu, sigma=sigma, verbose=TRUE)

  mrispatial26 <- makeMRIspatial(mask, nnei=26, sub=TRUE, bias=TRUE)
  tc.mcmcsubbias <- mritc.bayes(y, mrispatial$neighbors, mrispatial$blocks,
                                mrispatial$sub, mrispatial$subvox,
                                subbias=TRUE, mrispatial26$neighbors,
                                mrispatial26$blocks,mrispatial26$weineighbors,
                                mrispatial26$weights, mu=mu, sigma=sigma, verbose=TRUE)

  }
  #Example 2
  T1 <- readMRI(system.file("extdata/t1.rawb.gz", package="mritc"),
                c(91,109,91), format="rawb.gz")
  mask <-readMRI(system.file("extdata/mask.rawb.gz", package="mritc"),
                 c(91,109,91), format="rawb.gz")
  tc.icm <- mritc(T1, mask, method="ICM")
}
\keyword{classif}



