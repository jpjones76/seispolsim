#include "math.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* Declarations */
mxArray *xData;
mxArray *gData;
mxArray *hData;
double *Xt, *g, *W;
int M, L, M1, t, u, n;

/* Get matrix x */
Xt = mxGetPr(prhs[0]);
M = mxGetM(prhs[0]);
M1 = (M/2) - 1;

/* Get wavelet and scaling filters */
g = mxGetPr(prhs[1]);
L = mxGetM(prhs[1]);

/* Allocate memory and assign output pointer */
plhs[0] = mxCreateDoubleMatrix(1, M/2, mxREAL);
/* Note 01/02/2013 -- Ensure this is M, not M/2, for the modwpt routine */

/* Set pointers to the data space in our newly allocated memory */
W = mxGetPr(plhs[0]);

/* The algorithm */
for (t=0; t<=M1; t++) {
  u = (2*t)+1;
  W[t] = g[0]*Xt[u];
    for (n=1; n<=L-1; n++) {
	u = u-1;
	if (u<0) {u=M-1;}
	    W[t] += g[n]*Xt[u];
	}
  }
}
