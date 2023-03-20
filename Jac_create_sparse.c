/*==========================================================
 * Jac_crate.c - example in MATLAB External Interfaces
 *
 * Create jacobian matrix for RPNNode
 *
 * The calling syntax is:
 *
 *		[J1i,J1j,J1val,J2i,J2j,J2val] = Jac_create_sparse(trasp(Mmassi),trasp(Mmassj,trasp(Mmassval, PHIt_tk, DDi, DDj, DDval, PHIder)
 * This is a MEX file for MATLAB.
 * Copyright 2021 Gianluca Fabiani
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void Jac_create_sparse(double *J1i, double *J1j, double *J1val, double *J2i, double *J2j, double *J2val, double *Mmassi, double *Mmassj, double *Mmassval, double *PHIt_tk, double *DDi, double *DDj, double *DDval,double *PHIder, int N, int nnt, int nDD, int nMass){
int i, j, k, h, q, s, kis, kji;

/* nDD=nEq x nnt*nEq */
for (kis=0;kis<nDD;kis++){
    s=*(DDj+kis)-1;
    i=*(DDi+kis)-1;
    k=s%nnt;
    /*printf("kis=%d, s=%d, i=%d, k=%d \n",kis,s,i,k);*/
    for (h=0;h<N;h++){
        q=h+i*N;
        /*printf("kis=%d, s=%d, i=%d, h=%d, k=%d, q=%d, N=%d, nnt=%d \n",kis,s,i,h,k,q,N,nnt)*/
        /*J[q][s]=DD[i][s]*PHIt_tk[h][k];*/
        *(J1i + kis * N + h)=q+1;
        *(J1j + kis * N + h)=s+1;
        *(J1val + kis * N+h)=*(DDval+kis) * *(PHIt_tk+k*N+h);
    }
}
for (kji=0;kji<nMass;kji++){
    j=*(Mmassi+kji)-1;
    i=*(Mmassj+kji)-1;
    for (k=0;k<nnt;k++) {
        s=k+j*nnt; /* COLS */
        for(h=0;h<N;h++){
            q=h+N*i; /* ROWS */
            *(J2i+kji*N*nnt+h+N*k)=q+1;
            *(J2j+kji*N*nnt+h+N*k)=s+1;
            /*J[q][s]=PHIder[h][k]*Mmass[j][i];*/
            *(J2val+kji*N*nnt+h+N*k)=*(PHIder+k*N+h) * *(Mmassval+kji);
        }
    }
}

/*
for (j=0;j<nEq;j++) {
    for (k=0;k<nnt;k++) {
        s=k+j*nnt;
        for (i=0;i<nEq;i++) {
            for(h=0;h<N;h++){
                q=h+N*i;
                J[q][s]=DD[i][s]*PHIt_tk[h][k]+PHIder[h][k]*Mmass[j][i];
                 *(J+s+q*s2)=*(DD+s+i*s2) * *(PHIt_tk+k+h*nnt) + *(PHIder+k+h*nnt) * *(Mmass+i+j*nEq);
                *(J+s*s1+q)=*(DD+s*nEq+i) * *(PHIt_tk+k*N+h) + *(PHIder+k*N+h) * *(Mmass+i*nEq+j);
            }
        }
    }
}*/

}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
/* temp */
int i;

/* inputs */
double *Mmassi,*Mmassj,*Mmassval, *PHIt_tk, *DDi, *DDj, *DDval, *PHIder;

/* size of matrix */
int nnt, N, nDD, nMass, nJ1, nJ2;

/* output */
double *J1i,*J1j,*J1val,*J2i,*J2j,*J2val;      

if(nrhs != 8) {
    mexErrMsgIdAndTxt("MyToolbox:Jac_create:nrhs",
                      "Eight inputs required.");
}
if(nlhs!=6) {
    mexErrMsgIdAndTxt("MyToolbox:Jac_create:nlhs","Six output required.");
}

/* Get the dimensions of the input matrix  */
nMass=mxGetM(prhs[0]);
N = mxGetM(prhs[3]); /*rows*/
nnt = mxGetN(prhs[3]); /*columns*/
nDD = mxGetM(prhs[4]);


/* Create a pointer to the real data in the input matrix */
Mmassi = mxGetPr(prhs[0]); /* nEq x nEq */
Mmassj = mxGetPr(prhs[1]); /* nEq x nEq */
Mmassval = mxGetPr(prhs[2]); /* nEq x nEq */
PHIt_tk = mxGetPr(prhs[3]); /* N x nnt */
DDi = mxGetPr(prhs[4]); /* nEq x nnt*nEq */
DDj = mxGetPr(prhs[5]); /* nEq x nnt*nEq */
DDval = mxGetPr(prhs[6]); /* nEq x nnt*nEq */
PHIder = mxGetPr(prhs[7]); /* N x nnt */

/* how many non-zero? */


nJ1=N*nDD;
nJ2=N*nnt*nMass;

/* Set the output pointer to the output matrix. */
plhs[0]=mxCreateDoubleMatrix(nJ1,1, mxREAL);
plhs[1]=mxCreateDoubleMatrix(nJ1, 1, mxREAL);
plhs[2]=mxCreateDoubleMatrix(nJ1, 1, mxREAL);
plhs[3]=mxCreateDoubleMatrix(nJ2, 1, mxREAL);
plhs[4]=mxCreateDoubleMatrix(nJ2, 1, mxREAL);
plhs[5]=mxCreateDoubleMatrix(nJ2, 1, mxREAL);

/* get a pointer to the real data in the output matrix */
J1i = mxGetPr(plhs[0]);
J1j = mxGetPr(plhs[1]);
J1val = mxGetPr(plhs[2]);
J2i = mxGetPr(plhs[3]);
J2j = mxGetPr(plhs[4]);
J2val = mxGetPr(plhs[5]);

/* call the computational routine */
Jac_create_sparse(J1i, J1j, J1val, J2i, J2j, J2val, Mmassi, Mmassj, Mmassval, PHIt_tk, DDi, DDj, DDval, PHIder, N, nnt, nDD, nMass);

}/* end function */
