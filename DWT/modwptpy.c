#include "math.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Declarations */
  mxArray *xData;
  mxArray *gData;
  mxArray *hData;
  double *X, *g, *W;
  int N, L, j, t, n, k;
  
  /* Get matrix x */
  X = mxGetPr(prhs[0]);
  N = mxGetM(prhs[0]);
  
  /* Get wavelet and scaling filters */
  g = mxGetPr(prhs[1]);
  j = mxGetScalar(prhs[2]);
  L = mxGetM(prhs[1]);
  
  /* Create output */
  plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
  W = mxGetPr(plhs[0]);
  
  /* The algorithm */
  for (t = 0; t <= N-1; t++) {
    k = t;
    W[t] = g[0]*X[k];
    for (n = 1; n <= L-1; n++) {
      k -= pow(2,j-1);
      if (k<0) {k += N;}
      W[t] += g[n]*X[k];
    }
  }
}
