#include "math.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Declarations */
  double *W, *V, *g, *h, *X;
  int N, L, j, t, k, n;
  
  /* Get matrix x */
  W = mxGetPr(prhs[0]);
  N = mxGetNumberOfElements(prhs[0]);
  
  /* Get wavelet and scaling filters */
  h = mxGetPr(prhs[1]);
  j = mxGetScalar(prhs[2]);
  L = mxGetNumberOfElements(prhs[1]);
  
  /* Allocate memory and assign output pointer */
  plhs[0] = mxCreateDoubleMatrix(1, N, mxCOMPLEX);
  
  /* Set pointers to the data space in our newly allocated memory */
  X = mxGetPr(plhs[0]);
  
  /*printf("X is length %i, h is length %i.\n",N,L); 
    printf("j=%i\n",j);*/
  /* The algorithm */
  for (t=0; t<=N-1; t++) {
    k = t;
    X[t] = (h[0]*W[k]);
    for (n=1;n<=L-1;n++) {
      /* Correct, last rechecked Feb 2009 */
      k += pow(2,j-1);
      if (k>=N) {k = k%N;}
      X[t] += (h[n]*W[k]);
    }
  }
}
