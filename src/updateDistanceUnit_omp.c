#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* openMP is used*/
SEXP updateDistanceUnit(SEXP sbias, SEXP sneighbors, SEXP snvertex, SEXP snneigh, SEXP sweights, SEXP sweightsnew, SEXP sweineighbors, SEXP sweineighborsnew, SEXP sneiDiscrete)
{
	if (TYPEOF(sbias) != REALSXP)
		error("'bias' must be of type 'double'.");
	if (TYPEOF(sneighbors) != INTSXP)
		error("'neighbors' must be of type 'integer'.");
	if (TYPEOF(snvertex) != INTSXP)
		error("'nvertex' must be of type 'integer'.");
	if (TYPEOF(snneigh) != INTSXP)
		error("'nneigh' must be of type 'integer'.");
	if (TYPEOF(sweights) != REALSXP)
		error("'weights' must be of type 'double'.");
	if (TYPEOF(sweightsnew) != REALSXP)
		error("'weightsnew' must be of type 'double'.");
	if (TYPEOF(sweineighbors) != REALSXP)
		error("'weineighbors' must be of type 'double'.");
	if (TYPEOF(sweineighborsnew) != REALSXP)
		error("'weineighborsnew' must be of type 'double'.");
	if (TYPEOF(sneiDiscrete) != INTSXP)
		error("'neiDiscrete' must be of type 'integer'.");

    double *bias = REAL(sbias);
	int *neighbors = INTEGER(sneighbors);
    int nvertex = asInteger(snvertex);
    int nneigh = asInteger(snneigh);
    if (nneigh <= 0)
		error("The number of neighbors must be positive.");
	if (nvertex <= 0)
		error("The number of vertices must be positive.");
	double *weineighbors = REAL(sweineighbors);
	double *weineighborsnew = REAL(sweineighborsnew);
	double *weights = REAL(sweights);
	double *weightsnew = REAL(sweightsnew);
	int *neiDiscrete = INTEGER(sneiDiscrete);
  
	double *I1 = (double *) R_alloc(nvertex, sizeof(double));
	double *I2 = (double *) R_alloc(nvertex, sizeof(double));
	int i;

#pragma omp parallel for firstprivate(nneigh, nvertex, bias, neighbors, weights, \
									  weightsnew, weineighbors, weineighborsnew, \
									  I1, I2)
	for (i = 0; i < nvertex; i++) {
		int j;
		double sumB1 = 0.0, sumB2 = 0.0;
		for (j = 0; j < nneigh; j++) {
			int jj = i + j * nvertex;
			sumB1 = sumB1 + bias[neighbors[jj]]*weights[j]*neiDiscrete[jj];
			sumB2 = sumB2 + bias[neighbors[jj]]*weightsnew[j]*neiDiscrete[jj];
		}
		sumB1 = sumB1 + bias[i]*weineighbors[i];
		I1[i] = sumB1;
        sumB2 = sumB2 + bias[i]*weineighborsnew[i];
		I2[i] = sumB2;
	}

    SEXP qq;
	PROTECT(qq = allocVector(REALSXP, 2));
	for (i = 0; i < nvertex; i++) {
		REAL(qq)[0] = 	REAL(qq)[0] + I1[i]*bias[i];
		REAL(qq)[1] = 	REAL(qq)[1] + I2[i]*bias[i];
	}
	UNPROTECT(1);
    return qq;
}

