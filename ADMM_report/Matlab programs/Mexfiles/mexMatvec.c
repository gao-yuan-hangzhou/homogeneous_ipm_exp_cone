/***********************************************************************
* mexMatvec: compute A*y +z or AT*y +z.
* 
* Ay=mexMatvec(A,y,options,z)
*          
*  options = 0, compute A*y+z
*          = 1, compute (y'*A)' +z  
************************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

#if !defined(_WIN32)
#define dgemv dgemv_
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

/**********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    double   *A, *y;
     mwIndex  *irA, *jcA, *iry, *jcy; 
     double   *ytmp, *Ay, *z;
     int      isspA, isspy, isspz, options;

     int      m1, n1, m2, n2, j, jm1; 
     int      i, r, k, istart, iend, kstart, kend;
     double   tmp; 
      
     mwSize   M, N, LDA, INCX, INCY;
     double   ALPHA, BETA;
     char     *TRANS;
     ALPHA = 1.0;
     BETA  = 1.0;      
     INCX  = 1; 
     INCY  = 1; 

/* CHECK THE DIMENSIONS */

   if (mxIsCell(prhs[0]) | mxIsCell(prhs[1])) { 
       mexErrMsgTxt(" mexMatvec: A, x must be a double array"); }
   if (nrhs <2) {
       mexErrMsgTxt(" mexMatvec: must have at least 2 inputs"); }
   if (nlhs >1) { 
       mexErrMsgTxt("mexMatvec: requires 1 output argument"); }
  
   /***** assign pointers *****/
       A = mxGetPr(prhs[0]); 
       m1 = mxGetM(prhs[0]); 
       n1 = mxGetN(prhs[0]);
       isspA = mxIsSparse(prhs[0]);
       if (isspA) { irA = mxGetIr(prhs[0]);
	            jcA = mxGetJc(prhs[0]); }
       isspy = mxIsSparse(prhs[1]);
       m2 = mxGetM(prhs[1]); 
       n2 = mxGetN(prhs[1]);            
       if (n2 > 1) { 
 	      mexErrMsgTxt("mexMatvec: 2ND input must be a column vector"); }
       if (isspy) { 
          iry = mxGetIr(prhs[1]);
          jcy = mxGetJc(prhs[1]); 
          ytmp = mxGetPr(prhs[1]); 
          /***** copy ytmp to y *****/ 
          y = mxCalloc(m2,sizeof(double)); 
          kstart = jcy[0]; kend = jcy[1]; 
          for (k=kstart; k<kend; k++) { 
	      r = iry[k]; y[r] = ytmp[k]; }
       } else {
          y = mxGetPr(prhs[1]); 
       }              
       if (nrhs==2) { options = 0; } 
       else { options = (int) *mxGetPr(prhs[2]); 
       }        
       if (options == 0 & n1 != m2) {
          mexErrMsgTxt("mexMatvec: 1ST and 2ND input not compatible."); 
       } else if (options & m1 != m2) {
          mexErrMsgTxt("mexMatvec: 1ST and 2ND input not compatible."); 
       }

       /***** create return argument *****/
       if (options==0) { 
          plhs[0] = mxCreateDoubleMatrix(m1,1,mxREAL); 
       } else { 
          plhs[0] = mxCreateDoubleMatrix(n1,1,mxREAL); 
       }
       Ay = mxGetPr(plhs[0]); 
       if (nrhs==4) {
          isspz = mxIsSparse(prhs[3]);
          if (isspz) { mexErrMsgTxt("mexMatvec: 4TH input cannot be sparse."); }
          z = mxGetPr(prhs[3]);
          if (options==0) { 
             if (mxGetM(prhs[3])!=m1) {
                mexErrMsgTxt("mexMatvec: dimension of 4TH input is incorrect");
             }
             for (k=0; k<m1; k++) { Ay[k]=z[k]; }
          } else { 
             if (mxGetM(prhs[3])!=n1) {
                mexErrMsgTxt("mexMatvec: dimension of 4TH input is incorrect");
             }              
             for (k=0; k<n1; k++) { Ay[k]=z[k]; }
          }
       }
    /***** main body *****/
    LDA = m1; M = m1; N = n1;          
    if (options==0) {
       if (!isspA) {
          TRANS = "N"; 
          dgemv(TRANS, &M, &N, &ALPHA, A, &LDA, y, &INCX, &BETA, Ay, &INCY);
       } else {
          for (j=0; j<n1; j++){
             tmp = y[j];
             if (tmp != 0) {
                istart = jcA[j]; iend = jcA[j+1]; 
  	            for (i=istart; i<iend; i++) {
                    r = irA[i]; 
	                Ay[r] += tmp*A[i]; }
             }
          }
       }
    } else {
       if (!isspA) {
          TRANS = "T"; 
          dgemv(TRANS, &M, &N, &ALPHA, A, &LDA, y, &INCX, &BETA, Ay, &INCY);
       } else {
          for (j=0; j<n1; j++){
               istart = jcA[j]; iend = jcA[j+1]; 
               if (istart<iend) {
                  tmp = 0; 
  	              for (i=istart; i<iend; i++) {
                      r = irA[i]; 
	                  tmp += y[r]*A[i]; }
                  Ay[j] += tmp; 
               }
          }
       }
    }
 return;
}
/**********************************************************/

