#include "hessian.h"

extern "C" {

void printMatrix(const char* msg, double* arr, int rows, int cols)
{
    Rprintf("\n%s:\n", msg);
    for (int j=0; j < rows; ++j) {
        for (int i=0; i < cols; ++i)
            Rprintf("\t%f", arr[i*rows + j]);
        Rprintf("\n");
    }
}

inline void termwiseMult(double* v, double* u, double* prod, int len)
{
    for (int i=0; i < len; ++i) prod[i] = v[i]*u[i];
}

inline void rowSums(double* mat, int nrow, int ncol, double* ans, double* ones)
{
    int unity = 1;
    double d_one = 1.0, d_zero = 0.0;
    char NoTrans = 'N';

    dgemv_(&NoTrans, &nrow, &ncol, &d_one, mat, &nrow, ones, &unity, &d_zero,
           ans, &unity);
}

void computeNonGenHess(const size& sz, double* X, double* Y, double* wt,
                       double* probMat, double* baseProb, int ncores, double* H)
{
    int unity = 1;
    double d_one = 1.0, d_zero = 0.0;
    char Trans = 'T', NoTrans = 'N';

    // Setup pointers into Y to access blocks by choice
    double** Yc = NULL;
    if (sz.f) {    
        Yc = new double*[sz.K];
        for(int i=0; i < sz.K; ++i) Yc[i] = Y + sz.N*sz.f*i;
    }

    // Allocate scratch space for hessian and diagonal weight matrices
    double** scratchMat = new double*[ncores];
    double** scratchHess = new double*[ncores];
    int dim = ((sz.p >= sz.f) ? sz.p : sz.f); 
    for (int i=0; i < ncores; ++i) {
        scratchMat[i] = new double[sz.N*dim];
        scratchHess[i] = new double[dim*dim];
    }
   
    bool both = (sz.p > 0) && (sz.f > 0); 
    int nblocks = ((sz.p) ? (sz.K - 1) : 0) + ((sz.f) ? sz.K : 0);
    #ifdef SUPPORT_OPENMP
    #pragma omp parallel for num_threads(ncores) schedule(dynamic)
    #endif
    for (int ii=0; ii < nblocks; ++ii) {
        int threadID = 0;
        #ifdef SUPPORT_OPENMP
        // Grab workspace array for this thread
        threadID = omp_get_thread_num();
        #endif
        double* MtW = scratchMat[threadID];
        double* hess = scratchHess[threadID];
       
        // Local variables
        int m = ii;
        if (both) m = ((ii <= sz.K - 2) ? ii : ii - (sz.K - 1));
        double* mat_m = X;
        int cols_matm = sz.p;
        bool useBP_m = false;      // do we need to use base case probability ?
        if ((!both && sz.f) || (both && ii > sz.K - 2)) {
            if (both) useBP_m = (ii == sz.K - 1);
            else useBP_m = (ii == 0);
            mat_m = Yc[m]; 
            --m;
            cols_matm = sz.f;
        }
        // Loop over the upper triangular part of the Hessian
        for (int jj=ii; jj < nblocks; ++jj) {
            int n = jj;
            if (both) n = ((jj <= sz.K - 2) ? jj : jj - (sz.K - 1));
            double* mat_n = X;
            int cols_matn = sz.p;
            bool useBP_n = false;  // do we need to use base case probability ?
            if ((!both && sz.f) || (both && jj > sz.K - 2)) {
                if (both) useBP_n = (jj == sz.K - 1);
                else useBP_n = (jj == 0);
                mat_n = Yc[n];
                --n;
                cols_matn = sz.f;
            }
             
            // Note: We have X^t (and Y^t) & we need (M^t)*W*M (M = X or Y)
            // First compute: (mat_m^t)*W_mn
            int loc = sz.N * cols_matm;
            dcopy_(&loc, mat_m, &unity, MtW, &unity);
            double w, p_in, p_im;
            for (int i=0; i < sz.N; ++i) {
                // Compute weight
                p_im = (!useBP_m) ? probMat[m*sz.N + i] : baseProb[i];
                p_in = (!useBP_n) ? probMat[n*sz.N + i] : baseProb[i];
                w = (n == m) ? p_im*(1.0 - p_in) : -p_im * p_in;
                if (wt) w *= wt[i];  // multiply fweights
                for (int j=0; j < cols_matm; ++j) MtW[j + i*cols_matm] *= w;
            }
            // Now compute: hess = (mat_m^t * W_mn)*(mat_n^t)
            // dim(hess) = cols_matm * cols_matn
            dgemm_(&NoTrans, &Trans, &cols_matm, &cols_matn, &sz.N, &d_one,
                   MtW, &cols_matm, mat_n, &cols_matn, &d_zero, hess, &cols_matm); 
            // Copy block into full hessian's (ii, jj) block
            double* dest = H;
            if (both) {
                int offset = ( ((ii <= sz.K-2)? sz.p * ii : 
                                (sz.f * (ii - sz.K+ 1) + (sz.K - 1) * sz.p)) +
                  sz.nparams * ((jj <= sz.K-2)? sz.p * jj : 
                                (sz.f * (jj- sz.K + 1) + (sz.K - 1) * sz.p)) );
                dest += offset;
            } else
                dest += (ii*dim + sz.nparams*jj*dim);
            for (int i=0; i < cols_matn; ++i) {
                dcopy_(&cols_matm, hess + i*cols_matm, &unity,
                       dest + i*sz.nparams, &unity);
            }
            // Copy transpose of block to full hessian's (jj, ii) block
            if (jj == ii) continue; // Do nothing for diagonal blocks;
            dest = H;
            if (both) {
              dest += ( ((jj <= sz.K-2) ? sz.p * jj : 
                                   sz.f * (jj - sz.K + 1)+ (sz.K - 1) * sz.p) +
               sz.nparams*( (ii <= sz.K-2) ? sz.p * ii :
                                   sz.f * (ii-sz.K + 1)+ (sz.K - 1) * sz.p ) );
            } else
                dest += (jj * dim + sz.nparams * ii * dim);
            for (int c=0; c < cols_matn; ++c)
                for (int r=0; r < cols_matm; ++r)        
                    dest[c + r*sz.nparams] = hess[r + c*cols_matm];
        } // end inner for-loop
    } // end outer for-loop
    // Free allocated memory 
    if (sz.f) delete [] Yc;
    for (int i=0; i < ncores; ++i) {
        delete [] scratchMat[i];
        delete [] scratchHess[i];
    }
    delete [] scratchMat;
    delete [] scratchHess;
}

void computeGenHessCorner(const size& sz, double* Z, double* wt, double* probMat,
                          int ncores, double* H, double** scratchArrNK)
{
    int unity = 1;
    double d_one = 1.0, d_zero = 0.0, d_minusOne = -1.0;
    char Trans = 'T', NoTrans = 'N';

    int K = sz.K - 1; // nrows(Z) = (sz.K - 1) * sz.N
    int N = sz.N, d = sz.d;
  
    // Setup pointers to acces Z by blocks
    double**  Zc = new double*[sz.d];
    for(int i=0; i < sz.d; ++i) Zc[i] = Z + sz.N*K*i;

    double* ones_K = new double[K];
    for(int i=0; i < K; ++i) ones_K[i] = 1.0;

    double** scratchPZ = scratchArrNK;

    // Assemble hess and mat, column-by-column 
    double* hess = new double[sz.d * sz.d];
    double* mat = new double[sz.N * sz.d];
    ncores = (sz.d >= ncores) ? ncores : sz.d; // limit num_threads to sz.d
    #ifdef SUPPORT_OPENMP
    #pragma omp parallel for num_threads(ncores) 
    #endif 
    for (int ii=0; ii < sz.d; ++ii) {
        int threadID = 0;
        #ifdef SUPPORT_OPENMP
        // Grab workspace array for this thread
        threadID = omp_get_thread_num();
        #endif
        double* PZ = scratchPZ[threadID];

        // Compute: PZ_i[j] = prob[j] * Z_col_i[j] 
        // NOTE: dim(PZ) = N * K
        termwiseMult(probMat, Zc[ii], PZ, sz.N * K);

        // Compute col ii of matrix: mat = rowSums(PZ), dim(mat) = N x sz.d
        // NOTE: each col of 'mat' contains sum across choices
        rowSums(PZ, sz.N, K, mat + ii*sz.N, ones_K);

        if (wt) {
            // Multiply fweights: to each col of PZ multiply corresponding wt 
            for (int k=0; k < K; ++k)
                termwiseMult(PZ + k * sz.N, wt, PZ + k * sz.N, sz.N);
        }

        // Compute col ii of H: H[ , ii] = Z^T * PZ_ii 
        int nk = sz.N*K; 
        dgemv_(&Trans, &nk, &d, &d_one, Z, &nk, PZ, &unity, &d_zero,
               hess + ii*sz.d, &unity);
    }
    if (wt) {
        // Perform mat = sqrt(W) * mat
        // NOTE: sqrt is taken because we need: mat^T * W * mat
        #ifdef SUPPORT_OPENMP
        #pragma omp parallel for num_threads(ncores) 
        #endif 
        for (int j=0; j < sz.d; ++j) {
            for (int i=0; i < sz.N; ++i) 
                mat[i + j * sz.N] *= sqrt(wt[i]); 
        }
    }
    // Hessian block: hess = hess - mat^T * mat
    dgemm_(&Trans, &NoTrans, &d, &d, &N, &d_minusOne, mat, &N, mat, &N, &d_one,
           hess, &d);
    
    // Copy block into the full Hessian matrix 
    double* dest = H + sz.nparams*(sz.nparams - sz.d) + sz.nparams - sz.d;
    for (int i=0; i < sz.d; ++i)
        dcopy_(&d, hess + i*sz.d, &unity, dest + i *sz.nparams, &unity);
    
    delete [] Zc;
    delete [] ones_K;
    delete [] hess;
    delete [] mat;
  
}

void computeGenHess(const size& sz, double* X, double* Y, double* Z, double* wt,
                    double* probMat, double* baseProb, int ncores, double* H)
{
    int unity = 1;
    double d_one = 1.0, d_zero = 0.0;
    char NoTrans = 'N';

    int K = sz.K - 1; // nrows(Z) = (sz.K - 1) * sz.N
    int N = sz.N, d = sz.d;
    
    // Allocate memory needed by this function & computeGenHessCorner
    double** scratchArrNK = new double*[ncores];
    for (int i=0; i < ncores; ++i) 
        scratchArrNK[i] = new double[sz.N * K];
    
    // Compute the right-lower corner of Hessian (interaction among gen coeff)
    computeGenHessCorner(sz, Z, wt, probMat, ncores, H, scratchArrNK);
    
    if (sz.p == 0 && sz.f == 0) { //nothing more to be done
        for (int i=0; i < ncores; ++i) delete [] scratchArrNK[i];
        delete [] scratchArrNK;
        return;
    }
 
    // Setup pointers into Y & Z to access blocks by choice
    double** Yc = NULL; 
    if (sz.f) {    // if there are these types of vars
        Yc = new double*[sz.K];
        for(int i=0; i < sz.K; ++i) Yc[i] = Y + sz.N*sz.f*i;
    }
    double**  Zc = new double*[sz.d];
    for(int i=0; i < sz.d; ++i) Zc[i] = Z + sz.N*K*i;
     
    // Allocate scratch space for hessian and diagonal weight matrices
    double** scratchWt = new double*[ncores];
    double** scratchWZ = scratchArrNK;
    double** scratchHess = new double*[ncores];
    double** scratchMat = new double*[ncores];
    int dim = ((sz.p >= sz.f) ? sz.p : sz.f); 
    for (int i=0; i < ncores; ++i) {
        scratchWt[i] = new double[sz.N * K];
        scratchHess[i] = new double[dim * sz.d];
        scratchMat[i] = new double[sz.N * sz.d];
    }

    double* ones_K = new double[K];
    for(int i=0; i < K; ++i) ones_K[i] = 1.0;
    
    bool both = (sz.p > 0) && (sz.f > 0); 
    int nblocks = ((sz.p) ? (sz.K - 1) : 0) + ((sz.f) ? sz.K : 0);
    #ifdef SUPPORT_OPENMP
    #pragma omp parallel for num_threads(ncores) schedule(dynamic)
    #endif
    for (int ii=0; ii < nblocks; ++ii) {
        int threadID = 0;
        #ifdef SUPPORT_OPENMP
        // Grab workspace array for this thread
        threadID = omp_get_thread_num();
        #endif
        double* WZ = scratchWZ[threadID];
        double* W = scratchWt[threadID];
        double* hess = scratchHess[threadID];
        double* mat = scratchMat[threadID];

        // Local variables
        int n = ii;
        if (both) n = ((ii <= sz.K - 2) ? ii : ii - (sz.K - 1));
        double* mat_n = X;
        int cols_matn = sz.p;
        bool useBP_n = false;      // do we need to use base case probability ?
        if ((!both && sz.f) || (both && ii > sz.K - 2)) {
            if (both) useBP_n = (ii == sz.K - 1);
            else useBP_n = (ii == 0);
            mat_n = Yc[n];
            --n;
            cols_matn = sz.f;
        }

        // Compute the weight matrix: W
        // dim = N*K, entries for each individual contiguous
        for (int k=0; k < K; ++k) {
            double p_ik, p_in;
            for (int i=0; i < sz.N; ++i) {
                p_ik = probMat[k*sz.N + i];
                if (wt) p_ik *= wt[i];  // multiply fweights
                p_in = (!useBP_n) ? probMat[n*sz.N + i] : baseProb[i];
                W[i + k*sz.N] = (k == n) ? p_ik*(1.0 - p_in) : -p_ik * p_in;
            }
        }
        // Compute: mat_ia = sum_k {W_ik * Z_ika} 
        for (int jj=0; jj < sz.d; ++jj) {
            // Compute first: WZ = termwiseProd(array W, array Zc[j])
            termwiseMult(W, Zc[jj], WZ, sz.N * K);
            // Next: col j of mat = rowSums(WZ)
            rowSums(WZ, sz.N, K, mat + jj*sz.N, ones_K);
        }
        // Compute: hess = (mat_n)^T * mat, dim(hess) = cols_matn x sz.d
        // Note: We have direct access to (mat_n)^T
        dgemm_(&NoTrans, &NoTrans, &cols_matn, &d, &N, &d_one, mat_n,
               &cols_matn, mat, &N, &d_zero, hess, &cols_matn);

        // Copy this block into full Hessian's upper triangle
        double* dest = H + sz.nparams*(sz.nparams - sz.d);
        if (both) {
            dest += ( (ii <= sz.K - 2) ? ii * sz.p :
                                 sz.f * (ii - sz.K + 1) + (sz.K - 1) * sz.p );
        }
        else 
            dest += (ii * dim);  
        for (int i=0; i < sz.d; ++i) {
            dcopy_(&cols_matn, hess + i*cols_matn, &unity, dest + i*sz.nparams,
                   &unity); 
        }
        // Copy transpose of block into full Hessian's lower triangle
        dest = H + sz.nparams - sz.d;
        if (both) {
            dest += sz.nparams * ( (ii <= sz.K - 2) ?  ii * sz.p : 
                                 sz.f * (ii - sz.K + 1) + (sz.K - 1) * sz.p );
        } 
        else 
            dest += sz.nparams * ii * dim;
        for (int c=0; c < sz.d; ++c)
            for (int r=0; r < cols_matn; ++r)
                dest[c + r*sz.nparams] = hess[r + c*cols_matn];

    } // end parallel for-loop

    // Free allocated memory 
    if (sz.f) delete [] Yc;
    delete [] Zc;
    for (int i=0; i < ncores; ++i) {
        delete [] scratchArrNK[i];
        delete [] scratchWt[i];
        delete [] scratchHess[i];
        delete [] scratchMat[i];
    }
    delete [] scratchArrNK;
    delete [] scratchWt;
    delete [] scratchHess;
    delete [] scratchMat;
    delete [] ones_K;
}

SEXP computeHessianDotCall(SEXP Np, SEXP Kp, SEXP pp, SEXP fp, SEXP dp,
    SEXP Xm, SEXP Ym, SEXP Zm, SEXP Wv, SEXP probM, SEXP baseProbV,
    SEXP nprocs, SEXP hessM)
{
    int N  = INTEGER(Np)[0];
    int K  = INTEGER(Kp)[0];
    int p  = INTEGER(pp)[0];
    int f  = INTEGER(fp)[0];
    int d  = INTEGER(dp)[0];
    int np = p * (K - 1) + f * K + d; 

    double *X = NULL, *Y = NULL, *Z = NULL, *wt = NULL;
    if (!isNull(Xm)) X = REAL(Xm);
    if (!isNull(Ym)) Y = REAL(Ym);
    if (!isNull(Zm)) Z = REAL(Zm);
    if (!isNull(Wv)) wt = REAL(Wv);
    int hasWt = 0;
    if (wt) hasWt = 1;

    double *probMat = REAL(probM);
    double *baseProbVec = REAL(baseProbV);
    int ncores = INTEGER(nprocs)[0];
       
    //if (X) printMatrix("X matrix", X, N, p);

    double* H = REAL(hessM);
    computeHessian(&N, &K, &p, &f, &d, &hasWt, X, Y, Z, wt, probMat,
        baseProbVec, &ncores, H);

    //printMatrix("Hessian matrix", H, np, np);

    return R_NilValue;
}

void computeHessian(int* Np, int* Kp, int* pp, int* fp, int* dp, int* hasWt, 
            double* X, double* Y, double* Z, double* weight, double* probMat,
            double* baseProbVec, int* ncoresp, double* H)
{
    size sz(*Np, *Kp, *pp, *fp, *dp);
    double* wt = NULL;
    if (*hasWt) wt = weight;
    int ncores = *ncoresp;
    #ifndef SUPPORT_OPENMP
    if (ncores > 1) {
        ncores = 1;
        warning("OpenMP support not available for C++ compiler used with R.");
    }
    #endif    
    if (sz.p || sz.f)
        computeNonGenHess(sz, X, Y, wt, probMat, baseProbVec, ncores, H);
    if (sz.d)
        computeGenHess(sz, X, Y, Z, wt, probMat, baseProbVec, ncores, H);
}

// Extracts from an R style list, the element corresponding to tag
SEXP getListElement(SEXP list, const char *tag) {
    if (!isNewList(list)) error("Arg list isn't a newList.");

   SEXP elmnt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
   for (R_len_t i = 0; i < length(names); i++) {
       if (strcmp(tag, CHAR(STRING_ELT(names, i)))==0) {
           elmnt = VECTOR_ELT(list, i);
           break;
       }
   }
   return elmnt;
}


} // extern "C" 
