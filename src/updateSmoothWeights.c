#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP updateSmoothWeights(SEXP sgammaa, SEXP sgammab, SEXP snedges, SEXP sneighbors, SEXP snneigh, SEXP snvertex, SEXP sbias, SEXP sweights)
{
	if (TYPEOF(sgammaa) != REALSXP)
		error("'gammaa' must be of type 'double'.");
	if (TYPEOF(sgammab) != REALSXP)
		error("'gamma.b' must be of type 'double'.");
	if (TYPEOF(snedges) != INTSXP)
		error("'nedges' must be of type 'integer'.");
	if (TYPEOF(sneighbors) != INTSXP)
		error("'neighbors' must be of type 'integer'.");
	/*if (TYPEOF(sioneighbors) != INTSXP)
	  error("'ioneighbors' must be of type 'integer'.");*/
	if (TYPEOF(snneigh) != INTSXP)
		error("'nneigh' must be of type 'integer'.");
	if (TYPEOF(snvertex) != INTSXP)
		error("'nvertex' must be of type 'integer'.");
 	if (TYPEOF(sbias) != REALSXP)
		error("'bias' must be of type 'double'.");
 	if (TYPEOF(sweights) != REALSXP)
		error("'weights' must be of type 'double'.");

    
	double gammaa = asReal(sgammaa);
	double gammab = asReal(sgammab);
	int nedges = asInteger(snedges);
	int *neighbors = INTEGER(sneighbors);
	/*int *ioneighbors = INTEGER(sioneighbors);*/
    int nneigh = asInteger(snneigh);
    int nvertex = asInteger(snvertex);
    double *bias = REAL(sbias);
   	double *weights = REAL(sweights);

	SEXP val;
	PROTECT(val = allocVector(REALSXP, 1));
  
    GetRNGstate();
	
	int i, j;
	double sumD=0, shape, scale;
	
	for (i = 0; i < nvertex; i++) {
		for (j = 0; j < nneigh; j++) {
			double diff;
  
			diff = bias[i] - bias[neighbors[i + j * nvertex]-1];
			diff = diff * diff;
			/*sumD += (diff * weights[j]) * ioneighbors[i + j * nvertex];*/
			sumD += (diff * weights[j]);
		}
    }
	
	shape = gammaa + nedges / 2;
	scale = (2 * gammab) / (2 + gammab*sumD);
	REAL(val)[0] = rgamma(shape, scale);

    PutRNGstate();
	
	UNPROTECT(1);
    return val;
}

