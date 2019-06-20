#include<Rcpp.h>
#include<magma_v2.h>
#include<magma_lapack.h>
#include<magma_operators.h>

using namespace Rcpp;

RcppExport SEXP gpuQR_magma(SEXP X_)
{    
    // Input
    NumericMatrix X(X_);

    // Initialize magma and cublas
    magma_init();
 //   cublasInit();

    // Declare variables 
    int info, lwork, n_rows = X.nrow(), n_cols = X.ncol(), min_mn = n_rows;
    double tmp[1];
    double  work[1];
    NumericVector scale(min_mn);

    // Query workspace size
    magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]), &(work[0]), -1, &info); 
    Rcpp::Rcout << info << std::endl;

    lwork = work[0];
    NumericVector workN(lwork);

    // Run QR decomposition
    //             m       n      A        lda     tau           work       lwork   info
    // magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]), &(workN[0]), lwork, &info);
magma_dgeqrf2_gpu(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]),  &info);

    // Scale factor result    
    for(int ii = 1; ii < n_rows; ii++)
    {
        for(int jj = 0; jj < n_cols; jj++)
        {
            if(ii > jj) { X[ii + jj * n_rows] *= scale[jj]; }
        }
    }

    // Shutdown magma and cublas    
    magma_finalize();
   //  cublasShutdown();

    // Output  
    return wrap(X);
}


