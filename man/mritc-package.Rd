\name{mritc-package}
\alias{mritc-package}
\docType{package}
\title{MRI Tissue Classification Package}
\description{
  Use various methods to do MRI tissue classification.
}
\section{Introduction}{This package provides tools for MRI tissue
  classification using normal mixture models and hidden Markov normal
  mixture models (with the partial volume effect and intensity non-uniformity
  addressed) fitted by various methods.

  Magnetic resonance imaging (MRI) is used to
  identify the major tissues within a subject's brain. 
  Classification is usually based on a single image providing one
  measurement for each volume element, or voxel, in a discretization of
  the brain.  A simple model for MRI tissue classification views each voxel
  as homogeneous, belonging entirely to one of the three major tissue
  types (cerebrospinal fluid (CSF), gray matter (GM), and white matter
  (WM));  the intensity of voxels are thus normally
  distributed with means and variances depending on the tissue types
  of their voxels.
  The tissue types are not known and need to be identified from the
  image. The assumption that all tissue types are independent leads to
  a simple normal mixture model with parameters estimated by the EM
  algorithm and tissue types assigned using the Bayes classifier.
  
  Since nearby voxels tend to be of the same tissue type, a Markov
  random field model (a model from the Potts model family is used
  in this case) can be used to capture the spatial similarity of
  voxels by assigning homogeneity relationship among tissue types of
  neighboring voxels.
  Again, given the tissue types, the intensity of voxels are
  independently and normally distributed with means and
  variances depending on their tissue types.
  Furthermore, the Markov random field model defined on finite space
  is referred to as the hidden Markov model. Therefore the model
  combine the normal mixture part and the Potts model part is called
  the hidden Markov normal mixture model.
  This model can be fitted by the Iterated Conditional Mode algorithm,
  the Hidden Markov Random Field EM algorithm, or a Markov chain Monte
  Carlo approach. 

  A more realistic model than the one just described would take into
  account the fact that the volume elements are not homogeneous;
  while some may contain only one tissue type, others on the interface
  will contain two or possibly three different tissue types. This
  phenomenon is called the partial volume (PV) effect.
  One approach to address the PV effect is to introduce
  intermediate classes. Usually this is
  done by introducing two more classes: the combination of the CSF and
  the GM and the combination of the GM and the WM.
  Voxels containing WM and CSF are very rare and are ignored.
  This helps reduce confounding in estimation and a number of studies
  have used this approach.  Among these methods, the Gaussian partial volume
  hidden Markov random field models fitted by the modified EM algorithm
  appears to be more competitive in performance.
  A new approach to this problem is to construct a higher resolution
  image in which each voxel is divided into eight subvoxels.
  For each voxel
  the measured value is the sum of the unobserved measurements for the
  subvoxels. The subvoxels are in turn assumed to be homogeneous and
  follow the simpler model described above. 

  Intensity non-uniformity is an artifact that the signal intensity
  varies smoothly across an image. It is caused by combination
  and interaction of effects from the device, pulse sequence, and
  object. A commonly used approach to tackle it is to assume the 
  the measured signal is equal to true signal multiplied by bias field
  associated with the intensity non-uniformity plus some noise. 
  The bias field needs to be spatially smoothly varying and is modeled
  as either jointly normally distributed, or a linear combination of
  smooth spline or polynomial basis functions. Instead, we propose using
  a locally smoothed prior on the bias field.

  A Bayesian hierarchical model aiming at modeling the partial volume
  effect and intensity non-uniformity simultaneously was proposed.
  Instead of splitting the task into different steps, the framework
  harmoniously integrates several sub-models addressing different issues
  in the MRI classification, through specification of the likelihood function
  and prior distributions. This approach could provide more accurate
  tissue classification and also allow more effective estimation of the
  proportion of each voxel that belongs to each of the major tissue
  types.

  Besides brain image segmentation, the methods provided in this package can be used for 
  classification of other spatial data as well.
}
\section{Usage}{
  The function \code{\link[mritc]{readMRI}} and
  \code{\link[mritc]{writeMRI}} are I/O
  functions for MRI data. Right now, the "Analyze", "NIfTI", and raw
  byte (unsigned with 1 byte per element in the byte stream) gzip
  formats are supported.

  For each MR image, there has to be a corresponding array, \code{mask},
  with values 1 and 0. Voxels with value 1 are inside the brain and 0
  are outside.
  Tissue classification is conducted on voxels inside the brain.

  The functions \code{\link[mritc]{mritc.em}}, \code{\link[mritc]{mritc.icm}},
  \code{\link[mritc]{mritc.hmrfem}}, and
  \code{\link[mritc]{mritc.bayes}}
  are used
  to conduct the MRI tissue classification using the normal mixture
  model fitted by the EM algorithm, the hidden Markov normal
  mixture model at the voxel level fitted by the Iterated
  Conditional Mode algorithm, the Hidden Markov Random Field EM
  algorithm,
  or the Bayesian method (with or without the PV 
  or bias field correction).  The function
  \code{\link[mritc]{mritc.pvhmrfem}} is for classification using Gaussian
  partial volume hidden Markov random field models fitted by the
  modified EM algorithm.
  Different components of the normal mixture model correspond to
  different tissue types. The number of components is flexible,
  say using five components model to address the PV effect by
  \code{\link[mritc]{mritc.em}}, \code{\link[mritc]{mritc.icm}},
  \code{\link[mritc]{mritc.hmrfem}}, or \code{\link[mritc]{mritc.bayes}}.

  In order to use the previous functions, the parameters of the
  normal mixture model and the Potts model have to be specified.
  Some parameters can be obtained using the functions
  \code{\link[mritc]{initOtsu}} and \code{\link[mritc]{makeMRIspatial}}.
  There are default values for other parameters.

  The function \code{\link[mritc]{mritc}} integrates all methods together,
  provides a uniform platform with easier usage, and generates an object
  of class "mritc" , for which generic functions
  \code{\link[mritc]{print.mritc}}, \code{\link[mritc]{summary.mritc}},
  and \code{\link[mritc]{plot.mritc}} are provided. 
}
\section{Computation Issues}{
  To improve the speed, the table lookup method was used in various
  places; vectorized computation was used to take advantage of
  conditional independence. 
  Some computations are performed by \bold{C}
  code, and the \bold{OpenMP} is used to parallelize key loops in
  the \bold{C} code. Sparse matrix multiplication is adopted as well.
}
\keyword{package}
