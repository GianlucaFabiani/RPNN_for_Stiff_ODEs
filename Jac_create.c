/*==========================================================
 * Jac_crate.c - example in MATLAB External Interfaces
 *
 * Create jacobian matrix for RPNNode
 *
 * The calling syntax is:
 *
 *		J = Jac_create(Mmass, PHIt_tk, DD, PHIder)
 * This is a MEX file for MATLAB.
 * Copyright 2021 Gianluca Fabiani
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void Jac_create(double *J, double *Mmass, double *PHIt_tk, double *DD, double *PHIder, int nEq, int N, int nnt, int s1, int s2){
int i, j, k, h, q, s;


for (j=0;j<nEq;j++) {
    for (k=0;k<nnt;k++) {
        s=k+j*nnt; /* COLS */
        for (i=0;i<nEq;i++) {
            for(h=0;h<N;h++){
                q=h+N*i; /* ROWS */
                /*J[q][s]=DD[i][s]*PHIt_tk[h][k]+PHIder[h][k]*Mmass[j][i];*/
                /* *(J+s+q*s2)=*(DD+s+i*s2) * *(PHIt_tk+k+h*nnt) + *(PHIder+k+h*nnt) * *(Mmass+i+j*nEq); */
                *(J+s*s1+q)=*(DD+s*nEq+i) * *(PHIt_tk+k*N+h) + *(PHIder+k*N+h) * *(Mmass+i*nEq+j);
            }
        }
    }
}

}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
/* temp */
int i;

/* inputs */
double *Mmass, *PHIt_tk, *DD, *PHIder;

/* size of matrix */
int nEq, nnt, N, s1, s2;

/* output */
double *J;      

if(nrhs != 4) {
    mexErrMsgIdAndTxt("MyToolbox:Jac_create:nrhs",
                      "Four inputs required.");
}
if(nlhs!=1) {
    mexErrMsgIdAndTxt("MyToolbox:Jac_create:nlhs","One output required.");
}

/* Get the dimensions of the input matrix  */
nEq = mxGetN(prhs[0]);
N = mxGetM(prhs[1]); /*rows*/
nnt = mxGetN(prhs[1]); /*columns*/

s1=N*nEq;
s2=nnt*nEq;

/* Set the output pointer to the output matrix. */
plhs[0]=mxCreateDoubleMatrix(s1,s2, mxREAL);

/* Create a pointer to the real data in the input matrix */
Mmass = mxGetPr(prhs[0]); /* nEq x nEq */
PHIt_tk = mxGetPr(prhs[1]); /* N x nnt */
DD = mxGetPr(prhs[2]); /* nEq x nnt*nEq */
PHIder = mxGetPr(prhs[3]); /* N x nnt */

/* get a pointer to the real data in the output matrix */
J = mxGetPr(plhs[0]);

/* call the computational routine */
Jac_create(J, Mmass, PHIt_tk, DD, PHIder, nEq, N, nnt,s1,s2);
}

