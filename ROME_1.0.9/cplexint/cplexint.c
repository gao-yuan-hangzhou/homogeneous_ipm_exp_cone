/*
 *    CPLEXINT      MATLAB MEX INTERFACE FOR CPLEX, ver. 2.3
 *
 *    Copyright (C) 2001-2005  Mato Baotic
 *    Copyright (C) 2006       Michal Kvasnica
 *
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * file:        cplexint.c
 * 
 * Project:     MEX interface for CPLEX solver
 * 
 * Purpose:     Invoke the CPLEX callable lib.
 *
 * Authors:     M. Baotic, M. Kvasnica
 * 
 * Contact:     Mato Baotic
 *              Faculty of Electrical Engineering and Computing,
 *              Unska 3,
 *              HR-10000 Zagreb
 *              Croatia
 *              mato.baotic@fer.hr
 *
 *
 * History: date: yyyy.mm.mm | subject | (author)
 * ----------------------------------------------
 * 2008.02.14  Added socp for Quadractic constraint
 * 2006.10.03  Added time limit and cater for null vector input (Melvyn)
 * 2006.06.12  Interface addapted to CPLEX 10 (Michal)
 * 
 * 2005.09.25  Fixed memory leak when clearing CPXenv,
 *             Allow sparse matrix entries for quadratic constraints,
 *             OPTIONS.lic_rel introduced, (Mato)
 *             ver.2.3
 * 2005.09.21  Fix: mexErrMsgTxt() problems with R14 Matlab on Linux,
 *                  VARTYPE allocation, (Mato: thanks to Michal Kvasnica) 
 *             ver.2.2
 * 2004.03.06  Clean up of the code, (Mato)
 *             ver.2.1
 * 2004.05.17  Initial relase, (Mato)
 *             ver.2.0
 *
 * Notes        This file requires CPLEX (9.0 or higher)
 *              and MATLAB (6.0 or higher) to be compiled
 *                                                  
 *              CPLEX    http://www.ilog.com                     
 *              MATLAB   http://www.mathworks.com                
 *
 * (C) Sep 24 2005 by Mato Baotic
 * All rights reserved. 
 */
 
/* $Revision: 1.0 $ */

/*
  Matlab interface for CPLEX 9.0 solver for the following optimization problem
  
     min    0.5*x'*H*x + f'*x
      x
     s.t.:  A x {'<=' | '='} b
            x' * QC(i).Q * x + QC(i).L * x <= QC(i).r,  i=1,...,nQC
            LB <= x <= UB
            x(i) is {'C' = continuous | 'B' = binary | 'I' = integer |
                     'S' = semi-continuous | 'N' = semi-integer}, i=1,...,n
  
  The calling syntax is:
  [XOPT,FOPT,SOLSTAT,DETAILS] = cplexint(H, f, A, b, INDEQ, QC, LB, UB,...
                                         VARTYPE, PARAM, OPTIONS)
 */
 


/* MATLAB declarations. */
#include <stdlib.h>
#include <matrix.h>
#include "mex.h"

/* CPLEX declarations.  */
#include "cplex.h"


#define CPLEXINT_VERSION "2.3"
#define CPLEXINT_COPYRIGHT "Copyright (C) 2001-2005  Mato Baotic"

static CPXENVptr env = NULL;
static CPXFILEptr LogFile = NULL;

static int NUM_CALLS_CPLEXINT = 0;  /* number of calls to CPLEXINT before clearing
                                       environment CPXenv and releasing the license */
static int FIRST_CALL_CPLEXINT = 1; /* is this first call to CPLEXINT */


/* Problem type.
   We introduce our definition because getting the problem type
   from CPLEX env proved not to be the smartest idea.
   Maybe in the future CPLEX will handle this in a better way.
 */
enum {minLP, minQP, minMILP, minMIQP, minQCLP, minQCQP, minQCMILP, minQCMIQP};


/* MEX Input Arguments */
enum {H_IN_POS, F_IN_POS, A_IN_POS, B_IN_POS, INDEQ_IN_POS, QC_IN_POS,
      LB_IN_POS, UB_IN_POS, VARTYPE_IN_POS, PARAM_IN_POS, OPTIONS_IN_POS,
      MAX_NUM_IN_ARG};

#define MIN_NUM_IN_ARG      4

#define H_IN         prhs[H_IN_POS]
#define F_IN         prhs[F_IN_POS]
#define A_IN         prhs[A_IN_POS]
#define B_IN         prhs[B_IN_POS]
#define INDEQ_IN     prhs[INDEQ_IN_POS]
#define QC_IN        prhs[QC_IN_POS]
#define LB_IN        prhs[LB_IN_POS]
#define UB_IN        prhs[UB_IN_POS]
#define VARTYPE_IN   prhs[VARTYPE_IN_POS]
#define PARAM_IN     prhs[PARAM_IN_POS]
#define OPTIONS_IN   prhs[OPTIONS_IN_POS]


/* MEX Output Arguments */
enum {XMIN_OUT_POS, FMIN_OUT_POS, SOLSTAT_OUT_POS,
      DETAILS_OUT_POS, MAX_NUM_OUT_ARG};
      
/* If this is the only way to force people to check
   SOLSTAT before interpreting results then so be it! */
#define MIN_NUM_OUT_ARG     3

#define XMIN_OUT     plhs[XMIN_OUT_POS]
#define FMIN_OUT     plhs[FMIN_OUT_POS]
#define SOLSTAT_OUT  plhs[SOLSTAT_OUT_POS]
#define DETAILS_OUT  plhs[DETAILS_OUT_POS]


#define MAX_STR_LENGTH      1024


#if !defined(MAX)
#define MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

/* This hack is because Matlab R14 can crash on Linux due to the call to
   mexErrMsgTxt() */
#define TROUBLE_mexErrMsgTxt(A)    mexPrintf(A); mexPrintf("\n"); return
/* Here is the original version 
#define TROUBLE_mexErrMsgTxt(A)    mexErrMsgTxt(A)
*/

/* this is for TRYING to release the CPLEX license when pressing CTRL-C */
#define RELEASE_CPLEX_LIC   1

/*
    Copy and transform Matlab matrix IN in the representation needed for CPLEX.
    Return number of non zero elements (nnz) and CPLEX matrix description:
    matbeg, matcnt, matind, and matval.
    The arrays matbeg, matcnt, matind, and matval are accessed as follows.
    Suppose that CPLEX wants to access the entries in some column j. These
    are assumed to be given by the array entries:
        matval[matbeg[j]],.., matval[matbeg[j]+matcnt[j]-1]
    The corresponding row indices are:
        matind[matbeg[j]],.., matind[matbeg[j]+matcnt[j]-1]
    Entries in matind are not required to be in row order. Duplicate entries
    in matind within a single column are not allowed. The length of the arrays
    matbeg and matind should be of at least numcols. The length of arrays
    matind and matval should be of at least matbeg[numcols-1]+matcnt[numcols-1].
 */
int get_matrix(const mxArray *IN, int **outbeg, int **outcnt, int **outind,
               double **outval)
{
    int i, j;
    int gcount = 0;
    int pcount = 0;
    
    int m = mxGetM(IN);
    int n = mxGetN(IN);
    
    int    *matbeg = NULL;
    int    *matcnt = NULL;
    int    *matind = NULL;
    double *matval = NULL;
    
    double *in = mxGetPr(IN);


    matbeg = (int *)mxCalloc(n, sizeof(int));
    matcnt = (int *)mxCalloc(n, sizeof(int));

    /* Use different approaches for full and sparse matrix IN. */    
    if (!mxIsSparse(IN))
    {

        gcount = 0;
        for (i = 0; i < n; i++) {
            pcount = 0;
            for (j = 0; j < m; j++) {
                if (in[i * m + j] != 0) {
                    gcount++;
                    pcount++;
                }
                matbeg[i] = gcount - pcount;
                matcnt[i] = pcount;
            }
        }
        matind = (int *)mxCalloc(gcount, sizeof(int));
        matval = (double *)mxCalloc(gcount, sizeof(double));

        gcount = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                if (in[i * m + j] != 0) {
                    matind[gcount] = j;
                    matval[gcount] = in[i * m + j];
                    gcount++;
                }
            }
        }
    } else {
        /* For sparse matrix majority is already defined. */
        gcount = mxGetJc(IN)[n];
        matind = (int *)mxCalloc(gcount, sizeof(int));
        matval = (double *)mxCalloc(gcount, sizeof(double));
        
        for (i=0; i<n; i++){
            matbeg[i] = mxGetJc(IN)[i];
            matcnt[i] = mxGetJc(IN)[i+1] - mxGetJc(IN)[i];
        }
        for (i=0; i<gcount; i++){
            matind[i] = mxGetIr(IN)[i];
            matval[i] = in[i];
        }
    }
    *outbeg = matbeg;
    *outcnt = matcnt;
    *outind = matind;
    *outval = matval;
    return (gcount);
}

/*
    Copy and transform Matlab matrix IN in the vector representation needed for
    CPLEX. Return number of non zero elements (nnz) and CPLEX vector
    description: matval.
 */
int get_vector(const mxArray *IN, double **outval)
{
    int i;
    int gcount = 0;
    
    int m = mxGetM(IN);
    int n = mxGetN(IN);
    
    double *matval = NULL;
        
    double *in = mxGetPr(IN);

    matval = (double *)mxCalloc(m*n, sizeof(double));

    gcount = 0;
    for (i = 0; i < m*n; i++) {
        matval[i] = in[i];
        /* Handle infinity entries */
        /*
        if (matval[i]==mxGetInf()){
            matval[i]=CPX_INFBOUND;
        } else if (matval[i]==-mxGetInf()) {
            matval[i]=-CPX_INFBOUND;
        } 
        */
        if (in[i] != 0)
            gcount++;
    }
    *outval = matval;
    return (gcount);
}

/* 
    Display CPLEX error code message
 */
void dispCPLEXerror(CPXENVptr env, int status)
{
    char errmsg[MAX_STR_LENGTH];
    char *errstr;

    errstr = (char *)CPXgeterrorstring (env, status, errmsg);
    if ( errstr != NULL ) {
        mexPrintf("%s",errmsg);
    }else {
        mexPrintf("CPLEX Error %5d:  Unknown error code.\n", status);
    }   
}


/* 
   Here is the exit function, which gets run when the MEX-file is
   cleared and when the user exits MATLAB. The mexAtExit function
   should always be declared as static. 
 */
static void freelicence(void)
{
    int             status;
    extern CPXENVptr env;
    extern CPXFILEptr LogFile;
    extern int NUM_CALLS_CPLEXINT;  /* number of calls to CPLEXINT */
    extern int FIRST_CALL_CPLEXINT; /* is this first call to CPLEXINT */

    /* Close log file */
    if (LogFile != NULL){
        mexPrintf("LogFile is not NULL.\n");
    	status=CPXfclose(LogFile);

        if (status) {
            mexPrintf("Could not close log file cplexint_logfile.log.\n");
        } else {
            /* Just to be on the safe side we declare that the LogFile after
               closing is NULL. In this way we avoid possible error when trying
               to clear the same mex file more than once. */
            LogFile = NULL;
        }
    } else {
        /* mexPrintf("LogFile is NULL.\n"); */
    }
    

    /* Close CPLEX environment */
    FIRST_CALL_CPLEXINT = 1;
    NUM_CALLS_CPLEXINT = 0;
    if (env != NULL) {
        /* mexPrintf("env is not NULL.\n"); */
        status = CPXcloseCPLEX(&env);
        /*
           Note that CPXcloseCPLEX produces no output,
           so the only way to see the cause of the error is to use
           CPXgeterrorstring.  For other CPLEX routines, the errors will
           be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. 
         */
        if (status) {
            mexPrintf("Could not close CPLEX environment.\n");
            dispCPLEXerror(env, status);
        } else {
            /* Just to be on the safe side we declare that the environment after
               closing is NULL. In this way we avoid possible error when trying
               to clear the same mex file more than once. */
            env = NULL;
        }
    } else {
        /* mexPrintf("env is NULL.\n"); */
    }
    
}



/************************************
 *                                  *
 *   CPLEXINT solver MATLAB side    *
 *                                  *
 ************************************/
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    int             probtype = minLP;  /* default optimization problem type */
    
    char            errmsg[MAX_STR_LENGTH];     /* buffer for error messages */
    
     /* tmp variables */
    int             i, ii, jj, kk, mm;
    char           *tmpc;
    double         *tmpd;
    mxArray        *tmpArr;
    
    int             n_vars = 0;
    int             n_constr = 0;
    
    /* A_IN variables */
    int             A_nnz = 0;          /* number of non-zero elements */
    int            *A_matbeg = NULL;
    int            *A_matcnt = NULL;
    int            *A_matind = NULL;
    double         *A_matval = NULL;

    /* B_IN variables */
    int             b_nnz = 0;          /* number of non-zero elements */
    double         *b_matval = NULL;

    /* H_IN variables */
    int             H_nnz = 0;          /* number of non-zero elements */
    int            *H_matbeg = NULL;
    int            *H_matcnt = NULL;
    int            *H_matind = NULL;
    double         *H_matval = NULL;
    
    /* F_IN variables */
    int             f_nnz = 0;          /* number of non-zero elements */
    double         *f_matval = NULL;

    /* INDEQ_IN variables */
    int             sense_nnz = 0;      /* number of non-zero elements */
    char           *sense = NULL;       /* is constraint 'L' or 'E'    */

    /* Quadratic Constraint QC_IN variables */
    int             nQC = 0;            /* number of the QC */
    double         *QC_r = NULL;        /* right hand side of the QC */
    int            *QC_linnzcnt = NULL; /* linear part of the QC,
                                           number of nonzero elements */
    int           **QC_linind = NULL;   /* linear part of the QC,
                                           indices of nonzero elements */
    double        **QC_linval = NULL;   /* linear part of the QC,
                                           values of nonzero elements */
    int            *QC_quadnzcnt = NULL;/* quadratic part of the QC,
                                           number of nonzero elements */
    int           **QC_quadrow = NULL;  /* quadratic part of the QC,
                                           row indices of nonzero elements */
    int           **QC_quadcol = NULL;  /* quadratic part of the QC,
                                           column indices of nnz elements */
    double        **QC_quadval = NULL;  /* quadratic part of the QC,
                                           values of nonzero elements */


    /* LB_IN variables */
    int             LB_nnz = 0;          /* number of non-zero elements */
    double         *LB_matval = NULL;

    /* UB_IN variables */
    int             UB_nnz = 0;          /* number of non-zero elements */
    double         *UB_matval = NULL;
    
    /* VARTYPE_IN variables */
    int             vartype_nnC = 0;     /* number of non-continuous elements */
    char           *vartype = NULL;
    
    /* PARAM_IN variables */
    int             nintpar = 0;
    double         *intparcode = NULL;
    double         *intparvalue = NULL;
    int             ndoublepar = 0;
    double         *doubleparcode = NULL;
    double         *doubleparvalue = NULL;


    /* OPTIONS_IN variables */
    const char     *opt_fnames[] = {"verbose","save_prob","x0","probtype","lic_rel"};
    int             opt_nfields = (sizeof(opt_fnames)/sizeof(*opt_fnames));
    int             opt_verbose = 0;
    int             opt_save_probind = 0;
    char           *opt_save_prob = NULL;
    int             opt_logfile = 0;
    int             opt_nx0 = 0;
    int            *opt_x0i = NULL;
    double         *opt_x0 = NULL;
    int             opt_probtype = -1;  /* user can specify problem type,
                                           -1 internaly means: no specification
                                           by the user */
    int             opt_lic_rel = 1;    /* user can specify after how many calls will
                                           CPLEX environment be closed and license released */
	double          opt_tilimit = -1;  /* Time limit

    
    /* XMIN_OUT variables */
    mxArray        *XMIN = NULL;
    double         *xmin = NULL;
    
    /* FMIN_OUT variables */
    mxArray        *FMIN = NULL;
    double         *fmin = NULL;
    
    /* SOLSTAT_OUT variables */
    mxArray        *SOLSTAT = NULL;
    double         *solstat = NULL;


    /* DETAILS_OUT variables */
    mxArray        *DETAILS = NULL;
    const char     *details_fnames[] = {"statstring","solnmethod","solntype",
                                        "pfeasind","dfeasind","lpsolved",
                                        "dual","slack","qcslack","redcost"};
    int             details_nf = sizeof(details_fnames)/sizeof(*details_fnames);
    char           *statstring = NULL;
    int             solnmethod;
    int             solntype;
    int             pfeasind;
    int             dfeasind;
    int             lpsolved;
    double         *dual = NULL;
    double         *slack = NULL;
    double         *qcslack = NULL;
    double         *redcost = NULL;



    /* CPLEX variables */
    char            probname[] = "cplexint_problem\0";
    extern CPXENVptr env;
    extern CPXFILEptr LogFile;
    extern int NUM_CALLS_CPLEXINT;  /* number of calls to CPLEXINT before clearing
                                       environment CPXenv and releasing the license */
    extern int FIRST_CALL_CPLEXINT; /* is this first call to CPLEXINT */
    CPXLPptr        lp = NULL;
    int             status;
    int             cplex_probtype;

    int             errors = 1;     /* idicator of success */



    /* If there are no input nor output arguments display version number */
    if ((nrhs == 0) && (nlhs == 0)){
    	mexPrintf("CPLEXINT, Version %s.\n", CPLEXINT_VERSION);
    	mexPrintf("MEX interface for using CPLEX in Matlab.\n");
    	mexPrintf("%s.\n", CPLEXINT_COPYRIGHT);
    	return;
    }

    /* Check for proper number of arguments. */
    if (nrhs < MIN_NUM_IN_ARG) {
        sprintf(errmsg, "At least %d input arguments required.",
                MIN_NUM_IN_ARG);
        TROUBLE_mexErrMsgTxt(errmsg);
    } else if (nrhs > MAX_NUM_IN_ARG) {
        TROUBLE_mexErrMsgTxt("Too many input arguments.");
    } else if (nlhs < MIN_NUM_OUT_ARG) {
        mexPrintf("NOTE: This is a way of forcing people to notice SOLSTAT output.\n");
        mexPrintf("Always check SOLSTAT for correct inpterpretation of the results.\n");
        sprintf(errmsg, "At least %d output arguments required.",
                MIN_NUM_OUT_ARG);
        TROUBLE_mexErrMsgTxt(errmsg);
    } else if (nlhs > MAX_NUM_OUT_ARG) {
        TROUBLE_mexErrMsgTxt("Too many output arguments.");
    }

    /* Is somebody trying to solve an unconstrained problem. */
    if ((mxIsEmpty(A_IN)) || (mxIsEmpty(B_IN))){
        mexPrintf("If you are trying to solve an unconstrained problem one remedy is to\n");
        mexPrintf("artificially introduce some BIG constraints that will never be active.\n");
        TROUBLE_mexErrMsgTxt("CPLEX requires non-empty constraint matrices A and b.");
    }


    /* First get the number of variables and the number of constraints.
       For this purpose check size of a matrix A_IN (required argument). */
    if (   (!mxIsNumeric(A_IN))
        || (mxGetNumberOfDimensions(A_IN) > 2)
        || ((n_constr = mxGetM(A_IN)) < 1)
        || ((n_vars = mxGetN(A_IN)) < 1)
        || (mxIsComplex(A_IN))
        || (mxGetPr(A_IN) == NULL)
        ) {
        TROUBLE_mexErrMsgTxt("Matrix A must be a real valued (m x n) matrix, with m>=1, n>=1.");
    }
    A_nnz = get_matrix(A_IN, &A_matbeg, &A_matcnt, &A_matind, &A_matval);
    if (A_nnz == 0) {
        TROUBLE_mexErrMsgTxt("At least one element of constraint matrix A must be non-zero.");
    }
        

    if (   (!mxIsNumeric(B_IN))
        || (mxGetNumberOfDimensions(B_IN) > 2)
        || (mxGetM(B_IN) != n_constr)
        || (mxGetN(B_IN) != 1)
        || (mxIsComplex(B_IN))
        || (mxGetPr(B_IN) == NULL)
        ) {
        sprintf(errmsg, 
            "Vector b must be a real valued (%d x 1) vector.", n_constr);
        TROUBLE_mexErrMsgTxt(errmsg);
    }
    b_nnz = get_vector(B_IN, &b_matval);
    

    if (!mxIsEmpty(H_IN)){
        if (!mxIsNumeric(H_IN)
            || (mxGetNumberOfDimensions(H_IN) > 2)
            || (mxGetM(H_IN) != n_vars)
            || (mxGetN(H_IN) != n_vars)
            || (mxIsComplex(H_IN))
            || (mxGetPr(H_IN) == NULL)
            ) {
            sprintf(errmsg, 
                "If non-empty, objective H must be a real valued (%d x %d) matrix.",
                n_vars, n_vars);
            TROUBLE_mexErrMsgTxt(errmsg);
        }
        H_nnz = get_matrix(H_IN, &H_matbeg, &H_matcnt, &H_matind, &H_matval);
    }
    

    /* Is somebody trying to solve a feasibility problem. */
    if ((mxIsEmpty(F_IN)) && (H_nnz == 0)){
        mexPrintf("Note: If you are trying to solve a feasibility problem you could try\n");
        mexPrintf("setting all coefficients of f to zero. Note that in practice this\n");
        mexPrintf("may not the best approach (feasible space might be unbounded,...).\n");
        TROUBLE_mexErrMsgTxt("CPLEX requires non-empty objective vector f.");
    }

    if (   (!mxIsNumeric(F_IN))
        || (mxGetNumberOfDimensions(F_IN) > 2)
        || (mxGetM(F_IN) != n_vars)
        || (mxGetN(F_IN) != 1)
        || (mxIsComplex(F_IN))
        || (mxGetPr(F_IN) == NULL)
        ) {
        sprintf(errmsg,
            "Objective f must be a real valued (%d x 1) vector.",
            n_vars);
        TROUBLE_mexErrMsgTxt(errmsg);
    }
    f_nnz = get_vector(F_IN, &f_matval);


    /* Initially we assume that all constraints are of '<=' type */
    sense = mxCalloc(n_constr+1, sizeof(char));
    for (i = 0; i < n_constr; i++)
        sense[i] = 'L';
    sense[n_constr] = 0;

    /* Now check if there are some equality constraints */
    if ((nrhs > INDEQ_IN_POS) && (!mxIsEmpty(INDEQ_IN))) {
        if (   (!mxIsNumeric(INDEQ_IN))
            || (mxGetNumberOfDimensions(INDEQ_IN) > 2)
            || (MIN((mxGetM(INDEQ_IN)), (mxGetN(INDEQ_IN))) != 1)
            || (mxIsComplex(INDEQ_IN))
            ) {
            TROUBLE_mexErrMsgTxt("INDEQ must be a vector of indices of equality constraints.");
        } else {
            tmpd = mxGetPr(INDEQ_IN);
            for (i = 0; i < mxGetNumberOfElements(INDEQ_IN); i++){
                if ((tmpd[i]>=1) && (tmpd[i]<=n_constr)){
                    sense[(int)tmpd[i]-1] = 'E';
                } else {
                    TROUBLE_mexErrMsgTxt("Index in INDEQ points to a non-existing constraint.");
                }
            }
        }
    }
    sense_nnz = n_constr;


    if ((nrhs > QC_IN_POS) && (mxIsStruct(QC_IN))){
        
        if ((mxGetNumberOfDimensions(QC_IN) > 2)
            || (MIN(mxGetM(QC_IN),mxGetN(QC_IN)) > 1))
        {
            TROUBLE_mexErrMsgTxt("QC must be a STRUCTURE array.");
        }else{
            nQC = mxGetM(QC_IN)*mxGetN(QC_IN);
        }
        
        QC_r = (double *)mxCalloc(nQC, sizeof(double));
        QC_linnzcnt = (int *)mxCalloc(nQC, sizeof(int));
        QC_linind = (int **)mxCalloc(nQC, sizeof(int *));
        QC_linval = (double **)mxCalloc(nQC, sizeof(double *));
        QC_quadnzcnt = (int *)mxCalloc(nQC, sizeof(int));
        QC_quadrow = (int **)mxCalloc(nQC, sizeof(int *));
        QC_quadcol = (int **)mxCalloc(nQC, sizeof(int *));
        QC_quadval = (double **)mxCalloc(nQC, sizeof(double *));

        for (ii = 0; ii < nQC; ii++) {
            QC_quadnzcnt[ii] = 0;
            /* SOCP Type */
            if(((tmpArr = mxGetField(QC_IN, ii, "S"))  !=NULL))
            {
                /* Check dimension of input */
                if((mxGetM(tmpArr)==1) && (mxGetN(tmpArr)==2)&&
                ((tmpd = mxGetPr(tmpArr)) !=NULL))
                {
                    
                    QC_r[ii] = 0;
                    QC_linnzcnt[ii] = 0;
                    QC_linind[ii] = 0;
                    QC_linval[ii] = 0;
                    
                    QC_quadnzcnt[ii] = (int)tmpd[1]-(int)tmpd[0]+1;
                    QC_quadrow[ii] = (int *)mxCalloc(QC_quadnzcnt[ii], sizeof(int));
                    QC_quadcol[ii] = (int *)mxCalloc(QC_quadnzcnt[ii], sizeof(int));
                    QC_quadval[ii] = (double *)mxCalloc(QC_quadnzcnt[ii], sizeof(double));
            
                    QC_quadrow[ii][0] = (int)tmpd[0]-1;
                    QC_quadcol[ii][0] = (int)tmpd[0]-1;
                    QC_quadval[ii][0] = -1;
                    mm=1;
                    for (jj=(int)tmpd[0]; jj<(int)tmpd[1]; jj++){
                        QC_quadrow[ii][mm] = jj;
                        QC_quadcol[ii][mm] = jj;
                        QC_quadval[ii][mm] = 1;
                            mm++;
                    }
                    continue;
                
                }
                else{ 
                    sprintf(errmsg,
                    "S(%d) must have a 1 by 2 array", ii+1);
                    TROUBLE_mexErrMsgTxt(errmsg);
                }

            }
            else
            if (((tmpArr = mxGetField(QC_IN, ii, "Q"))  !=NULL)  &&
                (mxGetNumberOfDimensions(tmpArr) == 2) &&
                (mxGetM(tmpArr)==n_vars) &&
                (mxGetN(tmpArr)==n_vars) &&
                ((tmpd = mxGetPr(tmpArr)) !=NULL))
            {
            	if (!mxIsSparse(tmpArr)){
                    for (jj=0; jj<n_vars; jj++){
                        for (kk=0; kk<n_vars; kk++){
                            /* count nonzero entries */
                            if (tmpd[jj*n_vars+kk] !=0){
                                QC_quadnzcnt[ii] = QC_quadnzcnt[ii]+1;
                            }
                        }
                    }
                } else {
                    QC_quadnzcnt[ii] = mxGetJc(tmpArr)[n_vars];
                }
            }else{
                sprintf(errmsg,
                    "QC(%d) must have a %d by %d matrix Q as a field.",
                    ii+1,n_vars, n_vars);
                TROUBLE_mexErrMsgTxt(errmsg);
            }

            /* Now that we know the number of nonzero entries for this QC. */
            /* Let's allocate the space for it and fill it in.             */
            QC_quadrow[ii] = (int *)mxCalloc(QC_quadnzcnt[ii], sizeof(int));
            QC_quadcol[ii] = (int *)mxCalloc(QC_quadnzcnt[ii], sizeof(int));
            QC_quadval[ii] = (double *)mxCalloc(QC_quadnzcnt[ii], sizeof(double));
            
            if (!mxIsSparse(tmpArr)){
                mm = 0;
                for (jj=0; jj<n_vars; jj++){
                    for (kk=0; kk<n_vars; kk++){
                        /* MATLAB stores matrices columnwise !!! */
                        if (tmpd[jj*n_vars+kk] !=0){
                            QC_quadrow[ii][mm] = kk;
                            QC_quadcol[ii][mm] = jj;
                            QC_quadval[ii][mm] = tmpd[jj*n_vars+kk];
                            mm++;
                        }
                    }
                }
            } else {
            	mm = 0;
            	for (jj=0; jj<n_vars; jj++){
            	    for (kk=0; kk < mxGetJc(tmpArr)[jj+1] - mxGetJc(tmpArr)[jj]; kk++){
            	    	QC_quadrow[ii][mm] = mxGetIr(tmpArr)[mm];
            	    	QC_quadcol[ii][mm] = jj;
            	    	QC_quadval[ii][mm] = tmpd[mm];
            	    	mm++;
            	    }
            	}
            }
        

            QC_linnzcnt[ii] = 0;
            QC_linind[ii] = 0;
            QC_linval[ii] = 0;
            if((tmpArr = mxGetField(QC_IN, ii, "L"))  !=NULL)
            {           
                if ((mxGetNumberOfDimensions(tmpArr) == 2) &&
                    (mxGetM(tmpArr)==1) &&
                    (mxGetN(tmpArr)==n_vars) &&
                    ((tmpd = mxGetPr(tmpArr)) !=NULL) )
                {
                    for (jj=0; jj<n_vars; jj++){
                        /* count nonzero entries */
                        if (tmpd[jj] !=0){
                            QC_linnzcnt[ii] = QC_linnzcnt[ii] + 1;
                        }
                     }
                }else{

                    sprintf(errmsg,
                        "QC(%d) must have a 1 by %d vector L as a field.",
                        ii+1, n_vars);
                    TROUBLE_mexErrMsgTxt(errmsg);
                }
            

	            /* Now that we know the number of nonzero entries for this QC. */
		        /* Let's allocate the space for it and fill it in.             */
			    QC_linind[ii] = (int *)mxCalloc(QC_linnzcnt[ii], sizeof (int));
				QC_linval[ii] = (double *)mxCalloc(QC_linnzcnt[ii], sizeof (double));
	            mm = 0;
		        for (jj=0; jj<n_vars; jj++){
			        if (tmpd[jj] !=0){
				        QC_linind[ii][mm] = jj;
					    QC_linval[ii][mm] = tmpd[jj];
						mm++;
					}
				}
			}

            
            if (((tmpArr = mxGetField(QC_IN, ii, "r"))  !=NULL)  &&
                (mxGetNumberOfDimensions(tmpArr) == 2) &&
                (mxGetM(tmpArr)==1) &&
                (mxGetN(tmpArr)==1) &&
                ((tmpd = mxGetPr(tmpArr)) !=NULL) )
            {
                    QC_r[ii]=tmpd[0];
            }else{
                sprintf(errmsg,
                    "QC(%d) must have a 1 by 1 scalar r as a field.",
                    ii+1);
                TROUBLE_mexErrMsgTxt(errmsg);
            }
        } /* ii<nQC */
    } /* QC_IN_POS */ 

        

    if ((nrhs > LB_IN_POS) && (!mxIsEmpty(LB_IN))) {
        if (   (!mxIsNumeric(LB_IN))
            || (mxGetNumberOfDimensions(LB_IN) > 2)
            || (mxGetM(LB_IN) != n_vars)
            || (mxGetN(LB_IN) != 1)
            || (mxIsComplex(LB_IN))
            || (mxGetPr(LB_IN) == NULL)
            ) {
            sprintf(errmsg,
                "LB must be a real valued (%d x 1) column vector.",
                n_vars);
            TROUBLE_mexErrMsgTxt(errmsg);
        }
        LB_nnz = get_vector(LB_IN, &LB_matval);
    }
    if (LB_matval == NULL){
        LB_matval = mxCalloc(n_vars, sizeof(double));
        for (i=0; i<n_vars; i++){
            LB_matval[i] = -mxGetInf();
        }
        LB_nnz = n_vars;
    }
    

    if ((nrhs > UB_IN_POS) && (!mxIsEmpty(UB_IN))) {
        if (   (!mxIsNumeric(UB_IN))
            || (mxGetNumberOfDimensions(UB_IN) > 2)
            || (mxGetM(UB_IN) != n_vars)
            || (mxGetN(UB_IN) != 1)
            || (mxIsComplex(UB_IN))
            || (mxGetPr(UB_IN) == NULL)
            ) {
            sprintf(errmsg,
                "UB must be a real valued (%d x 1) column vector.",
                n_vars);
            TROUBLE_mexErrMsgTxt(errmsg);
        }
        UB_nnz = get_vector(UB_IN, &UB_matval);
    }
    if (UB_matval == NULL){
        UB_matval = mxCalloc(n_vars, sizeof(double));
        for (i=0; i<n_vars; i++){
            UB_matval[i] = mxGetInf();
        }
        UB_nnz = n_vars;
    }


    if ((nrhs > VARTYPE_IN_POS) && (!mxIsEmpty(VARTYPE_IN))) {
        if (   (!mxIsChar(VARTYPE_IN))
            || (mxGetNumberOfDimensions(VARTYPE_IN) > 2)
            || (mxGetM(VARTYPE_IN) != n_vars)
            || (mxGetN(VARTYPE_IN) != 1)
            ) {
            sprintf(errmsg,
                "VARTYPE must be a char valued (%d x 1) column vector.",
                n_vars);
            TROUBLE_mexErrMsgTxt(errmsg);
        } else {
            /* Allocate enough memory to hold the converted string. */
            vartype = mxCalloc(n_vars+1, sizeof (char));

            /* Copy the string data from string_array_ptr into buf. */
            if (mxGetString(VARTYPE_IN, vartype, n_vars+1) != 0) {
                TROUBLE_mexErrMsgTxt("Could not convert string data of VARTYPE.");
            }
        
            
            /* Checking if the input is made only of C B I S N. */
            vartype_nnC = 0;
            for (i = 0; i < n_vars; i++){
                if (vartype[i] != 'C'){
                    vartype_nnC++;
                }
                if (   (vartype[i] != 'C')
                    && (vartype[i] != 'B')
                    && (vartype[i] != 'I')
                    && (vartype[i] != 'S')
                    && (vartype[i] != 'N')){
                    TROUBLE_mexErrMsgTxt("VARTYPE must contain only C,B,I,S,N");
                }
            }
        }
    }
    if (vartype == NULL) {
        vartype = mxCalloc(n_vars + 1, sizeof(char));
        for (i = 0; i < n_vars; i++)
            vartype[i] = 'C';
        vartype[n_vars] = 0;
        vartype_nnC = 0;
    }


    if ((nrhs > PARAM_IN_POS) && (mxIsStruct(PARAM_IN))){
        if (   ((tmpArr = mxGetField(PARAM_IN, 0, "int"))  !=NULL)
            && (mxGetNumberOfDimensions(tmpArr) == 2)
            && (mxGetM(tmpArr)>=1)
            ) {
            if (mxGetN(tmpArr)==2){
                nintpar=mxGetM(tmpArr);
                
                intparcode = mxCalloc(nintpar, sizeof(double));
                intparvalue = mxCalloc(nintpar, sizeof(double));

                tmpd=mxGetPr(tmpArr);
                for (i=0; i<nintpar; i++){
                    intparcode[i]=tmpd[i];
                    intparvalue[i]=tmpd[i+nintpar];
                }
            }else{
                TROUBLE_mexErrMsgTxt("PARAM.int must be a matrix of size nintpar x 2.");
            }
        }
        if (   ((tmpArr = mxGetField(PARAM_IN, 0, "double")) !=NULL)
            && (mxGetNumberOfDimensions(tmpArr) == 2) 
            && (mxGetM(tmpArr)>=1)
            ) {
            if (mxGetN(tmpArr)==2){
                ndoublepar=mxGetM(tmpArr);

                doubleparcode = mxCalloc(ndoublepar, sizeof(double));
                doubleparvalue = mxCalloc(ndoublepar, sizeof(double));

                tmpd=mxGetPr(tmpArr);
                for (i=0; i<ndoublepar; i++){
                    doubleparcode[i]=tmpd[i];
                    doubleparvalue[i]=tmpd[i+ndoublepar];
                }
            }else{
                TROUBLE_mexErrMsgTxt("PARAM.double must be a matrix of size ndoublepar x 2.");
            }
        }
    }


    if (   (nrhs > OPTIONS_IN_POS) && (mxIsStruct(OPTIONS_IN)) ) {

        /* OPTIONS.verbose */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "verbose"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
        ){
            if (   (mxIsNumeric(tmpArr))
                && (mxGetNumberOfDimensions(tmpArr) == 2)
                && (mxGetM(tmpArr)==1)
                && (mxGetN(tmpArr)==1)
                && (!mxIsComplex(tmpArr))
            ) {
                tmpd = mxGetPr(tmpArr);
                if ((tmpd[0] != 0) && (tmpd[0] != 1) && (tmpd[0] != 2)){
                    TROUBLE_mexErrMsgTxt("OPTIONS.verbose must be 0, 1 or 2.");
                }
                opt_verbose = (int)tmpd[0];
            }else{
                TROUBLE_mexErrMsgTxt("OPTIONS.verbose must be 0, 1 or 2.");
            }
        }
        
        /* OPTIONS.timelimit */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "timelimit"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
        ){
            if (   (mxIsNumeric(tmpArr))
                && (mxGetNumberOfDimensions(tmpArr) == 2)
                && (mxGetM(tmpArr)==1)
                && (mxGetN(tmpArr)==1)
                && (!mxIsComplex(tmpArr))
            ) {
                tmpd = mxGetPr(tmpArr);
                if (((double)tmpd[0] <1 )){
                    TROUBLE_mexErrMsgTxt("OPTIONS.timelimit must be at least 1 sec.");
                }
                opt_tilimit = (double)tmpd[0];
            }else{
                TROUBLE_mexErrMsgTxt("OPTIONS.timelimit must be at least 1 sec.");
            }
        }

        /* OPTIONS.save_prob */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "save_prob"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
        ) {
            if (   (mxIsChar(tmpArr))
                && (mxGetNumberOfDimensions(tmpArr) == 2)
            ) {
                opt_save_probind = 1;
                i = mxGetNumberOfElements(tmpArr)+1;
                /* Allocate enough memory to hold the converted string. */
                opt_save_prob = mxCalloc(i, sizeof (char));

                /* Copy the string data */
                if (mxGetString(tmpArr, opt_save_prob, i) != 0) {
                    TROUBLE_mexErrMsgTxt("Could not convert string data.");
                }
            }else{
                TROUBLE_mexErrMsgTxt("OPTIONS.save_prob must be a string.");
            }
        }

        /* OPTIONS.logfile */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "logfile"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
        ){
            if (   (mxIsNumeric(tmpArr))
                && (mxGetNumberOfDimensions(tmpArr) == 2)
                && (mxGetM(tmpArr)==1)
                && (mxGetN(tmpArr)==1)
                && (!mxIsComplex(tmpArr))
            ) {
                tmpd = mxGetPr(tmpArr);
                if (   (tmpd[0] != 0)
                    && (tmpd[0] != 1)
                ){
                    TROUBLE_mexErrMsgTxt("OPTIONS.logfile must be 0 or 1.");
                }
                opt_logfile = (int)tmpd[0];
            }else{
                TROUBLE_mexErrMsgTxt("OPTIONS.logfile must be 0 or 1.");
            }
        }

        
        /* OPTIONS.x0 */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "x0"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
            && (mxGetNumberOfDimensions(tmpArr) == 2)
            && (!mxIsComplex(tmpArr))
        ) {
            if (mxGetN(tmpArr)==2){
                opt_nx0=mxGetM(tmpArr);
                
                opt_x0i = (int *)mxCalloc(opt_nx0, sizeof(int));
                opt_x0 = mxCalloc(opt_nx0, sizeof(double));

                tmpd=mxGetPr(tmpArr);
                for (i=0; i<opt_nx0; i++){
                    if (tmpd[i]<1 || tmpd[i]>n_vars) {
                        TROUBLE_mexErrMsgTxt("OPTIONS.x0 is indexing non-existing variable.");
                    }
                    opt_x0i[i]=(int)tmpd[i] - 1;
                    opt_x0[i]=tmpd[i+opt_nx0];
                }
            }else{
                TROUBLE_mexErrMsgTxt("OPTIONS.x0 must be a matrix of size nx0 x 2.");
            }
        }

        /* OPTIONS.probtype */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "probtype"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
        ){
            if (   (mxIsNumeric(tmpArr))
                && (mxGetNumberOfDimensions(tmpArr) == 2)
                && (mxGetM(tmpArr)==1)
                && (mxGetN(tmpArr)==1)
                && (!mxIsComplex(tmpArr))
            ) {
                tmpd = mxGetPr(tmpArr);
                if (   (tmpd[0] != -1)
                    && (tmpd[0] != minLP)
                    && (tmpd[0] != minQP)
                    && (tmpd[0] != minMILP)
                    && (tmpd[0] != minMIQP)
                    && (tmpd[0] != minQCLP)
                    && (tmpd[0] != minQCQP)
                    && (tmpd[0] != minQCMILP)
                    && (tmpd[0] != minQCMIQP)
                ){
                    TROUBLE_mexErrMsgTxt("Unknown problem type specified in OPTIONS.probtype.");
                }
                opt_probtype = (int)tmpd[0];
            }else{
                TROUBLE_mexErrMsgTxt("Unknown problem type specified in OPTIONS.probtype.");
            }
        }
        
        /* OPTIONS.lic_rel */
        if (   ((tmpArr = mxGetField(OPTIONS_IN, 0, "lic_rel"))  !=NULL)
            && (!mxIsEmpty(tmpArr))
        ){
            if (   (mxIsNumeric(tmpArr))
                && (mxGetNumberOfDimensions(tmpArr) == 2)
                && (mxGetM(tmpArr)==1)
                && (mxGetN(tmpArr)==1)
                && (!mxIsComplex(tmpArr))
            ) {
                tmpd = mxGetPr(tmpArr);
                if (   (tmpd[0] < 1)
                    || (tmpd[0] != (int)tmpd[0])
                ){
                    TROUBLE_mexErrMsgTxt("Wrong number of calls to CPLEXINT before license is released specified in OPTIONS.lic_rel.");
                }
                opt_lic_rel = (int)tmpd[0];
            }else{
                TROUBLE_mexErrMsgTxt("Unknown problem type specified in OPTIONS.probtype.");
            }
        }
        
    } /* OPTIONS_IN */




    /*
       Register safe exit.
     */
#ifdef RELEASE_CPLEX_LIC
    mexAtExit(freelicence);
#endif



    /***********************************
     *                                 *
     *  CPLEXINT solver CPLEX side     *
     *                                 *
     ***********************************/




    /* Initialize the CPLEX environment. */
    if (FIRST_CALL_CPLEXINT){
        env = CPXopenCPLEX(&status);
        FIRST_CALL_CPLEXINT = 0;
        NUM_CALLS_CPLEXINT = opt_lic_rel - 1;
    } else {
    	NUM_CALLS_CPLEXINT--;
    }

    /*
       If an error occurs, the status value indicates the reason for
       failure.  A call to CPXgeterrorstring will produce the text of
       the error message.  Note that CPXopenCPLEXdevelop produces no output,
       so the only way to see the cause of the error is to use
       CPXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  
     */
    if (env == NULL) {
        mexPrintf("Could not open CPLEX environment.\n");
        dispCPLEXerror(env, status);
        goto TERMINATE;
    }



    /* Create a log file. */
    if (opt_logfile){
    	/* 
    	   Open a LogFile to print out any CPLEX messages in there.
    	   We do this since Matlab does not execute printf commands
    	   in MEX files properly under Windows.
    	 */
        LogFile = CPXfopen("cplexint_logfile.log", "w");
        if (LogFile == NULL) {
            TROUBLE_mexErrMsgTxt("Could not open the log file cplexint_logfile.log.\n");
        }
        status = CPXsetlogfile(env, LogFile);
        if (status) {
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
    }

    /* Turn on output to the screen only if opt_verbose>=1. */
    if (opt_verbose>=1) {
        status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
        if (status) {
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
    }

    /* Create the problem. */
    lp = CPXcreateprob(env, &status, probname);

    /*
       A returned pointer of NULL may mean that not enough memory
       was available or there was some other problem.  In the case of 
       failure, an error message will have been written to the error 
       channel from inside CPLEX.  In this example, the setting of
       the parameter CPX_PARAM_SCRIND causes the error message to
       appear on stdout.  
     */
    if (lp == NULL) {
        mexPrintf("Failed to create LP.\n");
        dispCPLEXerror(env, status);
        goto TERMINATE;
    }


    /* Now copy the problem data into the lp. */
    status = CPXcopylp(env, lp, n_vars, n_constr, CPX_MIN, f_matval, b_matval,
               sense, A_matbeg, A_matcnt, A_matind, A_matval,
               LB_matval, UB_matval, NULL);

    if (status) {
        mexPrintf("Failed to copy problem data.\n");
        dispCPLEXerror(env, status);
        goto TERMINATE;
    }


    /* Copy quadratic objecive if one exists. */
    if (H_nnz > 0){
        status = CPXcopyquad (env, lp, H_matbeg, H_matcnt, H_matind, H_matval);
        if (status) {
            mexPrintf("Failed to copy quadratic objective matrix H.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
    }


	 /* Set maximum computational time
     */
     if (opt_tilimit>=1) {
            status = CPXsetdblparam(env, CPX_PARAM_TILIM,opt_tilimit);
        if (status) {
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
    }

    
    /* Now copy the vartype array if it has some non-continuous components . */
    if (vartype_nnC > 0){
        status = CPXcopyctype(env, lp, vartype);
        if (status) {
            mexPrintf("Failed to copy VARTYPE.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        } 
    }

    /* Set initial solution. */
    if ((opt_nx0 > 0) && (vartype_nnC > 0)) {
        status = CPXcopymipstart(env, lp, opt_nx0, opt_x0i, opt_x0);
        if (status) {
            mexPrintf("Failed to set initial solution.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
        /* modified by M. Kvasnica. according to CPLEX10 manual,
         * CPX_PARAM_ADVIND must be set to 1
         */
        status = CPXsetintparam(env, CPX_PARAM_ADVIND, 1);
    }


    /* Add quadratic constraints to the problem. */
    for (i = 0; i<nQC; i++){
        status = CPXaddqconstr (env, lp, QC_linnzcnt[i], QC_quadnzcnt[i],
                                QC_r[i], 'L', QC_linind[i], QC_linval[i],
                                QC_quadrow[i], QC_quadcol[i], QC_quadval[i],
                                NULL);
        if (status) {
            mexPrintf("Failed to copy quadratic constraint.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
    }


    /* If verbose>=2 display problem on the screen. */
    if (opt_verbose>=2){
        FILE *fp;
        char buf[80];

        mexPrintf("\nVERBOSITY IS ON!\n");

        /* Save problem to the verbosity file for display on the screen. */
        status = CPXwriteprob(env, lp, "cplexint_verbose.lp", NULL);
        if (status) {
            mexPrintf("Failed to save the problem to the verbosity file.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
        mexPrintf("\nThe problem is stored in the file cplexint_verbose.lp\n");
        if ((fp = fopen("cplexint_verbose.lp", "r")) == NULL){
            mexPrintf("Failed to open verbosity file.");
            goto TERMINATE;
        }
          
        while (fgets(buf, sizeof(buf), fp) != NULL)
            mexPrintf("%s", buf);
        fclose(fp);

        /* Initial guess */
        if ((opt_nx0) && (vartype_nnC > 0)) {
            mexPrintf("\nStarting from:\n");
            for (i = 0; i < opt_nx0; i++) {
                mexPrintf("x%d = %15.10e\n", opt_x0i[i]+1, opt_x0[i]);
            }
        }
    }


    /* Save problem to the file. */
    if (opt_save_probind){
        status = CPXwriteprob(env, lp, opt_save_prob, NULL);
        if (status) {
            mexPrintf("Failed to save problem to the file.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }
    }


    /* The size of the problem should be obtained by asking CPLEX what
       the actual size is, rather than using what was passed to CPXcopylp.
       cur_numcols stores the current number of columns.  */
    if (opt_verbose>=2){
        mexPrintf("\nFrom input arguments we see that the problem has:\n");
        mexPrintf("%d variables,\n",n_vars);
        mexPrintf("%d constraints,\n",n_constr);
        mexPrintf("%d non-zero elements in quadratic objective,\n",H_nnz);
        mexPrintf("%d quadratic constraints,\n",nQC);
        mexPrintf("%d non-continuous variables,\n",vartype_nnC);
    }


    if (opt_probtype != -1){
        switch (opt_probtype){
            case (minLP):
                status = CPXchgprobtype(env,lp,CPXPROB_LP);
                break;
            case (minQP):
                status = CPXchgprobtype(env,lp,CPXPROB_QP);
                break;
            case (minMILP):
                status = CPXchgprobtype(env,lp,CPXPROB_MILP);
                break;
            case (minMIQP):
                status = CPXchgprobtype(env,lp,CPXPROB_MIQP);
                break;
            case (minQCLP):
                status = CPXchgprobtype(env,lp,CPXPROB_QCP);
                break;
            case (minQCQP):
                status = CPXchgprobtype(env,lp,CPXPROB_QCP);
                break;
            case (minQCMILP):
                status = CPXchgprobtype(env,lp,CPXPROB_MIQCP);
                break;
            case (minQCMIQP):
                status = CPXchgprobtype(env,lp,CPXPROB_MIQCP);
                break;
        }
        if (status) {
            mexPrintf("Failed to change problem type.\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        } else {
            if (opt_verbose>=2)
                mexPrintf("\nProblem type successfully changed to %d\n",
                    opt_probtype);
        }
    }
    
    
    n_vars = CPXgetnumcols(env,lp);
    n_constr = CPXgetnumrows(env,lp);
    if (H_nnz > 0){
        H_nnz = CPXgetnumqpnz(env,lp);
    }
    if (nQC > 0){
        nQC = CPXgetnumqconstrs(env,lp);
    }
    if (vartype_nnC > 0){
        vartype_nnC = CPXgetnumbin(env,lp) + CPXgetnumint(env,lp) + 
                      CPXgetnumsemicont(env,lp) + CPXgetnumsemiint(env,lp);
    }
    
    if (opt_verbose>=2){
        mexPrintf("\nCPLEX says that the problem has:\n");
        mexPrintf("%d variables,\n",n_vars);
        mexPrintf("%d constraints,\n",n_constr);
        mexPrintf("%d non-zero elements in quadratic objective,\n",H_nnz);
        mexPrintf("%d quadratic constraints,\n",nQC);
        mexPrintf("%d non-continuous variables,\n",vartype_nnC);
    }


    /*
       Create matrices for the three main return arguments.
     */
    XMIN = mxCreateDoubleMatrix(n_vars, 1, mxREAL);
    xmin = mxGetPr(XMIN);


    FMIN = mxCreateDoubleMatrix(1, 1, mxREAL);
    fmin = mxGetPr(FMIN);


    SOLSTAT = mxCreateDoubleMatrix(1, 1, mxREAL);
    solstat = mxGetPr(SOLSTAT);




    /* 
       Unless user specifies (OVERRIDES) the problem type through 
       OPTIONS.probtype we will decide based on the passed inputs.
     */
    if (opt_probtype != -1){
        probtype = opt_probtype;
    } else {
        if ((H_nnz == 0) && (vartype_nnC == 0) && (nQC == 0)){
            probtype = minLP;
        } else if ((H_nnz != 0) && (vartype_nnC == 0) && (nQC == 0)){
            probtype = minQP;
        } else if ((H_nnz == 0) && (vartype_nnC != 0) && (nQC == 0)){
            probtype = minMILP;
        } else if ((H_nnz != 0) && (vartype_nnC != 0) && (nQC == 0)){
            probtype = minMIQP;
        } else if ((H_nnz == 0) && (vartype_nnC == 0) && (nQC != 0)){
            probtype = minQCLP;
        } else if ((H_nnz != 0) && (vartype_nnC == 0) && (nQC != 0)){
            probtype = minQCQP;
        } else if ((H_nnz == 0) && (vartype_nnC != 0) && (nQC != 0)){
            probtype = minQCMILP;
        } else if ((H_nnz != 0) && (vartype_nnC != 0) && (nQC != 0)){
            probtype = minQCMIQP;
        }
    }



    
    
    /* Optimize the problem and obtain solution. */
    cplex_probtype = CPXgetprobtype(env, lp);
    if (opt_verbose>=2)
    {
        mexPrintf("\nCPLEX says that we have the following problem type: %d.\n",
            cplex_probtype);
        mexPrintf("For interpretation of the returned code check function\n");
        mexPrintf("CPXgetprobtype in Callable library (CPLEX Ref. manual).\n");
        mexPrintf("\nCPLEXINT says we have the following problem type: %d.\n",
            probtype);
        mexPrintf("For interpretation of returned code check cplexint.m.\n");
    }
    


    /* Change default values of CPLEX parameters. */
    for (i = 0; i < nintpar; i++){
        status = CPXsetintparam(env, (int)intparcode[i], (int)intparvalue[i]);
        if (status) {
            mexPrintf("Failed to set int parameter(s).\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }

        /* Display parameter change on the screen. */
        if (opt_verbose>=2){
            int             ivalue;
            CPXgetintparam(env, (int)intparcode[i], &ivalue);
            mexPrintf("CPLEX parameter %d is set to the value %d.\n", 
                (int)intparcode[i], ivalue);
        }
    }
    for (i = 0; i < ndoublepar; i++){
        status = CPXsetdblparam(env,
                                (int)doubleparcode[i], 
                                (double)doubleparvalue[i]);
        if (status) {
            mexPrintf("Failed to set double parameter(s).\n");
            dispCPLEXerror(env, status);
            goto TERMINATE;
        }

        /* Display parameter change on the screen. */
        if (opt_verbose>=2){
            double           dvalue;
            CPXgetdblparam(env, (int)doubleparcode[i], &dvalue);
            mexPrintf("CPLEX parameter %d is set to the value %15.10e.\n",
                (int)doubleparcode[i], dvalue);
        }
    }

    /* If there is at least one quadratic constraint we have to use CPXbaropt.
       This might overwrite user's request passed through the argument PARAM. */
    if (probtype == minLP){
        status = CPXlpopt(env, lp);
    } else if (probtype == minQP){
        status = CPXqpopt(env, lp);
    } else if (probtype == minMILP){
        status = CPXmipopt(env, lp);
    } else if (probtype == minMIQP){
        status = CPXmipopt(env, lp);
    } else if (probtype == minQCLP){
        status = CPXbaropt(env, lp);
    } else if (probtype == minQCQP){
        status = CPXbaropt(env, lp);
    } else if (probtype == minQCMILP){
        status = CPXmipopt(env, lp);
    } else if (probtype == minQCMIQP){
        status = CPXmipopt(env, lp);
    }else{
        mexPrintf("Unknown problem type.\n");
        goto TERMINATE;
    }

    /* If in the following functions we fail to get certain parameter,
       we don't want to interupt MEX file with an error in Matlab. */
    errors = 0;

    if ((status) && (opt_verbose>=1)) {
        mexPrintf("\n");
        dispCPLEXerror(env, status);
        *solstat = (double)CPXgetstat(env, lp);
        mexPrintf("\nStatus of the solution: %d\n", (int)(*solstat));
        tmpc = CPXgetstatstring(env, (int)(*solstat), errmsg);
        if (tmpc != NULL){
            mexPrintf("%s\n",errmsg);
        }
        status = CPXgetsubstat(env, lp);
        mexPrintf("\nStatus of the last subproblem: %d\n", status);
        tmpc = CPXgetstatstring(env, status, errmsg);
        if (tmpc != NULL){
            mexPrintf("%s\n",errmsg);
        }
        mexPrintf("\nFailed to optimize with CPLEXINT!!!\n\n");

        /*
          CPXwriteprob(env, lp, "cplexint_error.lp", NULL);
          mexPrintf("\nProblem is saved to the file cplexint_error.lp");
        */
    } else {
        /* Get status of the solution. */
        *solstat = (double)CPXgetstat(env, lp);
    }
    

    /* Get objective value and optimizer */
    if ((probtype == minMILP) || (probtype == minMIQP) ||
        (probtype == minQCMILP) || (probtype == minQCMIQP))
    {
        status = CPXgetmipx(env, lp, xmin, 0, n_vars - 1);
        if ((status) && (opt_verbose>=1)) {
            mexPrintf("Failed to get optimizer x.\n");
            dispCPLEXerror(env, status);
        }
        status = CPXgetmipobjval(env, lp, fmin);
        if ((status) && (opt_verbose>=1)) {
            mexPrintf("No MIP objective value available.\n");
            dispCPLEXerror(env, status);
        }
    } else {
        status = CPXgetx(env, lp, xmin, 0, n_vars - 1);
        if ((status) && (opt_verbose>=1)) {
            mexPrintf("Failed to get optimal x.\n");
            dispCPLEXerror(env, status);
        }

        status = CPXgetobjval(env, lp, fmin);
        if ((status) && (opt_verbose>=1)) {
            mexPrintf("No objective value available.\n");
            dispCPLEXerror(env, status);
        }
    }
    /* Write the XMIN, FMIN and SOLSTAT to the screen. */
    if (opt_verbose>=2)
    {
        mexPrintf("\nSolution:\n");
        for (i=0; i<n_vars; i++)
            mexPrintf("x%d = %f\n", i+1, xmin[i]);
        mexPrintf("\nObjective value  = %f\n", *fmin);
        mexPrintf("\nSolution status = %d\n", (int)solstat[0]);
    }


    if (nlhs <= DETAILS_OUT_POS)
        goto TERMINATE;
        
    /* Gather DETAILS */
    DETAILS = mxCreateStructMatrix(1, 1, details_nf, details_fnames);

    /* get solstat string */
    tmpc = CPXgetstatstring(env, (int)(*solstat), errmsg);
    if (tmpc != NULL){
        tmpArr = mxCreateString(errmsg);
        mxSetField(DETAILS, 0, "statstring", tmpArr);
    } else {
    	if (opt_verbose>=1)
    		mexWarnMsgTxt("Failed to get DETAILS (statstring).");
    }

    /* get solnmethod, solntype, pfeasible, dfeasible */
    if ((probtype == minLP) || (probtype == minQP) || (probtype == minQCLP) || 
        (probtype == minQCQP)){
        status = CPXsolninfo(env, lp, &solnmethod, &solntype, &pfeasind,
                             &dfeasind);
        if (status){
       	    if (opt_verbose>=1){
            	mexWarnMsgTxt("Failed to get DETAILS (solnmethod, solntype, pfeasible, dfeasible).");
            	dispCPLEXerror(env, status);
            }
        }else{
            tmpArr = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxGetPr(tmpArr)[0] = (double)solnmethod;
            mxSetField(DETAILS, 0, "solnmethod", tmpArr);
    
            tmpArr = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxGetPr(tmpArr)[0] = (double)solntype;
            mxSetField(DETAILS, 0, "solntype", tmpArr);

            tmpArr = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxGetPr(tmpArr)[0] = (double)pfeasind;
            mxSetField(DETAILS, 0, "pfeasind", tmpArr);

            tmpArr = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxGetPr(tmpArr)[0] = (double)dfeasind;
            mxSetField(DETAILS, 0, "dfeasind", tmpArr);
        }
    }
    
    /* get dual variables */
    if ((probtype == minLP) || (probtype == minQP)){
        tmpArr = mxCreateDoubleMatrix(n_constr, 1, mxREAL);
        dual = mxGetPr(tmpArr);
        status = CPXgetpi(env, lp, dual, 0, n_constr-1);
        if (status){
            if (opt_verbose>=1){
            	mexWarnMsgTxt("Failed to get DETAILS (dual).");
            	dispCPLEXerror(env, status);
       	    }
            if (tmpArr != NULL)
                mxDestroyArray(tmpArr);
        } else {
            mxSetField(DETAILS, 0, "dual", tmpArr);
        }
    }
    
    /* get slack variables. slacks exist for all problems */
    tmpArr = mxCreateDoubleMatrix(n_constr, 1, mxREAL);
    slack = mxGetPr(tmpArr);
    if ((probtype == minMILP) || (probtype == minMIQP) ||
        (probtype == minQCMILP) || (probtype == minQCMIQP)){
        status = CPXgetmipslack(env, lp, slack, 0, n_constr-1);
    } else {
        status = CPXgetslack(env, lp, slack, 0, n_constr-1);
    }
    if (status){
    	if (opt_verbose>=1){
            mexWarnMsgTxt("Failed to get DETAILS (slack).");
            dispCPLEXerror(env, status);
        }
        if (tmpArr != NULL)
            mxDestroyArray(tmpArr);
    } else {
        mxSetField(DETAILS, 0, "slack", tmpArr);
    }

    /* get quadratic constraints slack */
    if ((probtype == minQCLP) || (probtype == minQCQP) ||
        (probtype == minQCMILP) || (probtype == minQCMIQP)){
        if (nQC > 0){
            tmpArr = mxCreateDoubleMatrix(nQC, 1, mxREAL);
            qcslack = mxGetPr(tmpArr);
            if ((probtype == minQCMILP) || (probtype == minQCMIQP)){
                status = CPXgetmipqconstrslack(env, lp, qcslack, 0, nQC-1);
            } else {
                status = CPXgetqconstrslack(env, lp, qcslack, 0, nQC-1);
            }
            if (status){
    		if (opt_verbose>=1){
    		    mexWarnMsgTxt("Failed to get DETAILS (qcslack).");
                    dispCPLEXerror(env, status);
                }
                if (tmpArr != NULL)
                    mxDestroyArray(tmpArr);
            } else {
                mxSetField(DETAILS, 0, "qcslack", tmpArr);
            }
        }
    }
    
    /* get number of LP's being solved in MIP */
    if ((probtype == minMILP) || (probtype == minMIQP) ||
        (probtype == minQCMILP) || (probtype == minQCMIQP)){
        if (vartype_nnC > 0){
            lpsolved = CPXgetmipitcnt(env, lp);
            tmpArr = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxGetPr(tmpArr)[0] = (double)lpsolved;
            mxSetField(DETAILS, 0, "lpsolved", tmpArr);
        }
    }
    
    /* get reduced cost */
    if ((probtype == minLP) || (probtype == minQP)){
        tmpArr = mxCreateDoubleMatrix(n_vars, 1, mxREAL);
        redcost = mxGetPr(tmpArr);
        status = CPXgetdj(env, lp, redcost, 0, n_vars-1);
        if (status) {
            if (opt_verbose>=1){
                mexWarnMsgTxt("Failed to get DETAILS (redcost).");
                dispCPLEXerror(env, status);
            }
            if (tmpArr != NULL)
                mxDestroyArray(tmpArr);
        } else {
            mxSetField(DETAILS, 0, "redcost", tmpArr);
        }
    }


    
    
  TERMINATE:
  
    /* Close log file */
    if ((opt_logfile) && (LogFile != NULL)){
    	status = CPXfclose(LogFile);

        if (status) {
            mexPrintf("Could not close log file cplexint_logfile.log.\n");
        } else {
            /* Just to be on the safe side we declare that the LogFile after
               closing is NULL. In this way we avoid possible error when trying
               to clear the same mex file afterwards. */
            LogFile = NULL;
        }
    }

    
    /*
       Free up the problem as allocated by CPXcreateprob, if necessary.
     */
    if (lp != NULL) {
        status = CPXfreeprob(env, &lp);
        if (status) {
            if (opt_verbose>=1){
                mexPrintf("CPXfreeprob failed.\n");
                dispCPLEXerror(env, status);
            }
        }
    }

    /* Free up the CPLEX environment, if necessary. */
    if (NUM_CALLS_CPLEXINT <= 0){
        FIRST_CALL_CPLEXINT = 1; /* prepare for the next call to CPLEXINT */
        NUM_CALLS_CPLEXINT = 0;  /* prepare for the next call to CPLEXINT */

        if (env != NULL) {
            status = CPXcloseCPLEX(&env);
            
            /*
               Note that CPXcloseCPLEX produces no output,
               so the only way to see the cause of the error is to use
               CPXgeterrorstring.  For other CPLEX routines, the errors will
               be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. 
             */
            if (status) {
                if (opt_verbose>=1){
                    mexPrintf("Could not close CPLEX environment.\n");
                    dispCPLEXerror(env, status);
                }
            } else {
                /* Just to be on the safe side we declare that the environment after
                   closing is NULL. In this way we avoid possible error when trying
                   to clear the same mex file afterwards. */
                env = NULL;
            }
        }
    }

    if (!errors){
        /* Pass computation to the real outputs and clear not used memory */
        if (nlhs <= DETAILS_OUT_POS){
            if (DETAILS != NULL){
                mxDestroyArray(DETAILS);
            }
        } else {
            DETAILS_OUT = DETAILS;
        }
        SOLSTAT_OUT = SOLSTAT;
        FMIN_OUT = FMIN;
        XMIN_OUT = XMIN;
    } else {
        if (DETAILS != NULL)
            mxDestroyArray(DETAILS);
        if (SOLSTAT != NULL)
            mxDestroyArray(SOLSTAT);
        if (FMIN != NULL)
            mxDestroyArray(FMIN);
        if (XMIN != NULL)
            mxDestroyArray(XMIN);
    }
    

   
    /* Free allocated memory. */
    if (A_matbeg != NULL)
        mxFree(A_matbeg);
    if (A_matcnt != NULL)
        mxFree(A_matcnt);
    if (A_matind != NULL)
        mxFree(A_matind);
    if (A_matval != NULL)
        mxFree(A_matval);

    if (b_matval != NULL)
        mxFree(b_matval);

    if (H_matbeg != NULL)
        mxFree(H_matbeg);
    if (H_matcnt != NULL)
        mxFree(H_matcnt);
    if (H_matind != NULL)
        mxFree(H_matind);
    if (H_matval != NULL)
        mxFree(H_matval);

    if (f_matval != NULL)
        mxFree(f_matval);

    if (sense != NULL)
        mxFree(sense);

    if (nQC > 0){
        int ii;
        for (ii=0; ii<nQC; ii++){
            if (QC_linind[ii] != NULL)
                mxFree(QC_linind[ii]);
            if (QC_linval[ii] != NULL)
                mxFree(QC_linval[ii]);
            if (QC_quadrow[ii] != NULL)
                mxFree(QC_quadrow[ii]);
            if (QC_quadcol[ii] != NULL)
                mxFree(QC_quadcol[ii]);
            if (QC_quadval[ii] != NULL)
                mxFree(QC_quadval[ii]);
        }
    }
    if (QC_r != NULL)
        mxFree(QC_r);
    if (QC_linnzcnt != NULL)
        mxFree(QC_linnzcnt);
    if (QC_linind != NULL)
        mxFree(QC_linind);
    if (QC_linval != NULL)
        mxFree(QC_linval);
    if (QC_quadnzcnt != NULL)
        mxFree(QC_quadnzcnt);        
    if (QC_quadrow != NULL)
        mxFree(QC_quadrow);
    if (QC_quadcol != NULL)
        mxFree(QC_quadcol);
    if (QC_quadval != NULL)
        mxFree(QC_quadval); 

    if (LB_matval != NULL)
        mxFree(LB_matval);

    if (UB_matval != NULL)
        mxFree(UB_matval);

    if (intparcode != NULL)
        mxFree(intparcode);

    if (intparvalue != NULL)
        mxFree(intparvalue);

    if (doubleparcode != NULL)
        mxFree(doubleparcode);

    if (doubleparvalue != NULL)
        mxFree(doubleparvalue);

    if (vartype != NULL)
        mxFree(vartype);

    if (opt_save_prob != NULL)
        mxFree(opt_save_prob);

    if (opt_x0i != NULL)
        mxFree(opt_x0i);

    if (opt_x0 != NULL)
        mxFree(opt_x0);


    if (errors) {
        TROUBLE_mexErrMsgTxt("There were errors.");
    }
    return;
}

