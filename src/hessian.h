/*
* Parallelized implementation of Hessian computation for multinomial logit
* loglikelihood functions.
*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

extern "C" {

struct size {
  int N;       // Number of observations
  int K;       // Number of choices
  int p;       // Number of individual specifc variables
  int f;       // Number of choice sp vars with choice sp coeff
  int d;       // Number of choice sp vars with generic coeff
  int nparams; // Number of parameters to be estimated
 
  size(int NN, int KK, int pp, int ff, int dd):
    N(NN), K(KK), p(pp), f(ff), d(dd), nparams((K - 1)*p + K*f + d) {}
};

/*
* Prints to R a matrix, given as an array in col-major order
*/
void printMatrix(const char* msg, double* arr, int rows, int cols);

/*
* Computes: prod[i] = v[i] * u[i]  (0 <= i < len)
* Uses BLAS level 2 call: dsbmv
* Treats 'v' as a diagonal matrix (doesn't require extra storage).
*/
void termwiseMult(double* v, double* u, double* prod, int len);

/*
* Sums a matrix along its its rows. Uses level 2 BLAS call DGEMV.
*
* Input:
*   mat - array having matrix in col major order of size: nrow x ncol
*   ones - array containing 1.0, dim = ncol, 
*   ans - array of length = nrow
*
* Output:
*   At exit - ans = mat * ones (this effects a rowsum) 
*             No other array is changed nor any extra storage needed.
*/
void rowSums(double* mat, int nrow, int ncol, double* ans, double* ones);

/*
* Computes the non-generic coefficients block of the Hessian matrix.
* See docs of computeHessian() for explanation of args 
*/
void computeNonGenHess(const size& sz, double* X, double* Y, double* wt,
                     double* probMat, double* baseProb, int ncores, double* H);

/*
* Computes interaction amongst only generic coefficients (lower right corner)
* See docs of computeHessian() for explanation of args except:
*   scratchArrNK - Workspace array having 'ncores' sub-arrays of len = len(Z)
*/
void computeGenHessCorner(const size& sz, double* Z, double* wt, double* probMat,
                          int ncores, double* H, double** scratchArrNK);

/*
* Computes the generic coefficients block of the Hessian matrix.
* See docs of computeHessian() for explanation of args 
*/
void computeGenHess(const size& sz, double* X, double* Y, double* Z, double* wt, 
                    double* probMat, double* baseProb, int ncores, double* H);

/*
* ======= .Call Interface function, will be called from R. =======
* Input:
*    Np, Kp, pp, fp, dp: see 'size' in calling code for info 
*    Xm     : matrix with individual specific data 
*    Ym     : matrix with indv & choice varying data modeled with ind sp coeff 
*    Zm     : matrix with indv & choice varying data modeled with generic coeff 
*    Wv     : vector of frequency weights 
*    probM  : matrix of probabilities ( N * K-1) 
*  baseProbV: vector of probabilities of base choice (length N) 
*    nprocs : int - number of procs to run on 
*    hessM  : vector - pre-allocated to appropriate size 
*
* Output:
*    hessM is set at exit with the Hessian in col-major format
*
* Notes:
*   Expected format of matrices: See docs of computeHessin()
*/
SEXP computeHessianDotCall(SEXP Np, SEXP Kp, SEXP pp, SEXP fp, SEXP dp,
    SEXP Xm, SEXP Ym, SEXP Zm, SEXP Wv, SEXP probM, SEXP baseProbV,
    SEXP nprocs, SEXP hessM);

/*
* ======= .C Interface function, will be called from R. =======
* Input:
*    Np   : pointer to the integer containing the number of observations
*    Kp   : pointer to the integer containing the number of choices
*    pp   : pointer to the integer with # individual specific variables
*    fp   : pointer to the integer with # choice sp variables choice coeff
*    dp   : pointer to the integer with # choice specific variables gen coeff
*    hasWt: pointer to the integer which is positive if weight is NOT null
*    X    : pointer to the invidual-specific design matrix. 
*    Y    : pointer to the choice-specific vars & coeffs design matrix. 
*    Z    : pointer to the choice-specific vars & generic coeff design matrix. 
*    wt   : pointer to the vector of frequency weights of length N 
*    probMat     : pointer to the Probability matrix (P_ik, k in [1, K-1])   
*    baseProbVec : pointer to the base alt Probability vector (P_i0)   
*    ncoresp     : pointer to integer with # of processors allowed 
*    H           : pointer to the Hessian matrix initialzed to zeros. 
*
* Output:
*    On exit, Hessian matrix is stored in H.
*
* Notes:
*   Expected format of matrices:
*     X - entries for a single individual are contiguous. dim = p*N
*     Y - dim = f*N*K. Organised as K blocks of length f*N each. Entries for 
*           a single individual are contiguous.
*     Z - dim = N*(K-1)*d. d blocks of length N*(K-1). Entries for
*           a single choice are contiguous. Data for base choice absent.
*     wt - length = N, all entries must be positive
*     probMat - dim = N*(K-1). Entries for each choice are contiguous.
*     H - dim = nparams*nparams. Symmetric matrix.
*         First (K-1)*p rows are for individual specific vars. 
*         Next f rows correpond to the base choice and choice specific vars
*           with choice specific coefficients. 
*         Next (K-1)*f rows are for the rest of the choice specific vars
*           with choice sp coeffs.
*         Last d rows are for ch sp vars with generic coeffs.
*         When any of these groups isn't in the model, then Hessian rows
*           corresponding to it will not appear.
*   All matrices are in COLUMN MAJOR ORDERING.   
*/       
void computeHessian(int* Np, int* Kp, int* pp, int* fp, int* dp, int* hasWt, 
            double* X, double* Y, double* Z, double* wt, double* probMat, 
            double* baseProbVec, int* ncoresp, double* H);
/*
* Helper function from the manual Rextensions.pdf
*
* Extracts from an R style list, the element corresponding to tag
*/
SEXP getListElement(SEXP list, const char *tag);

} // extern "C" 
