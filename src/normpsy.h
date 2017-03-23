#include <R_ext/RS.h>

void F77_SUB(backtransformation)(double * mu,
				 double * VC0, 
				 double * VC1, 
				 int * maxmes,
				 double * spl,
				 int * nbzitr,
				 double * zitr0,
				 int * nsim,
				 int * methInteg,
				 double * Ymarg);

