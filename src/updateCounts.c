#include <R.h>
#include <Rinternals.h>

/* Destructively updates the counts  sZsum */
SEXP updateCounts(SEXP sZsum, SEXP sZ)
{
	if (TYPEOF(sZsum) != INTSXP)
		error("'Zsum' must be of type 'integer'.");
	if (TYPEOF(sZ) != INTSXP)
		error("'Z' must be of type 'integer'.");
	
	int n = LENGTH(sZ);
	if (n <= 0)
		error("The length of 'Z' must be positive.");
	if (n != LENGTH(sZsum))
		error("The length of 'Z' and 'Zsum' do not match.");
	int *Zsum = INTEGER(sZsum);
    int *Z = INTEGER(sZ);
    int i;
	
	
    for (i = 0; i < n; i++)
		Zsum[i] += Z[i];
	
    return sZsum;
}
