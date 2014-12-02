#include "mex.h"
#include <string.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {

	if(nlhs != 1)	mexErrMsgTxt("Check output number");
	if(nrhs != 6)	mexErrMsgTxt("Check input number");
	if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || mxIsInt32(prhs[4]) || mxIsInt32(prhs[5]) ) {
		mexErrMsgTxt("Check input type.");
	}

	// INPUT
	double *edgeInd1, *edgeInd2;
	double *fE1, *fE2;
	int nV1, nV2, nE1, nE2, nMatch, dFeat;
	edgeInd1 = mxGetPr(prhs[0]);
	edgeInd2 = mxGetPr(prhs[1]);
	fE1 = mxGetPr(prhs[2]);
	fE2 = mxGetPr(prhs[3]);
	nV1 = mxGetScalar(prhs[4]);
	nV2 = mxGetScalar(prhs[5]);
	nE1 = mxGetM(prhs[0]);
	nE2 = mxGetM(prhs[1]);
	nMatch = nV1*nV2;
	dFeat = mxGetM(prhs[2]);

	// OUTPUT
	double *affinityMatrix;
	plhs[0] = mxCreateDoubleMatrix(nMatch, nMatch, mxREAL);
	affinityMatrix = mxGetPr(plhs[0]);
	memset(affinityMatrix, 0, nMatch*nMatch*sizeof(double));	// TODO: make sparse 

	int i, j, a, b;
	for(int iE1 = 0; iE1 < nE1; iE1++) {
		for(int iE2 = 0; iE2 < nE2; iE2++) {
			i = edgeInd1[iE1]-1;
			j = edgeInd1[iE1+nE1]-1;
			a = edgeInd2[iE2]-1;
			b = edgeInd2[iE2+nE2]-1;
			for(int iFeat = 0; iFeat < dFeat; iFeat++) {
				affinityMatrix[i+a*nV1 + (j+b*nV1)*nMatch]
				+= fE1[iFeat+iE1*dFeat] * fE2[iFeat+iE2*dFeat];
			}
		}
	}

}