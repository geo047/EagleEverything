#include "Rcpp.h"
#include<magma_v2.h>


#undef CUBLASAPI
#ifdef __CUDACC__
#define CUBLASAPI __host__
#else
#define CUBLASAPI
#endif

#include "cublas_api.h"

#if defined(__cplusplus)
extern "C" {
#endif

/* CUBLAS data types */
#define cublasStatus cublasStatus_t

cublasStatus CUBLASWINAPI cublasInit (void);
cublasStatus CUBLASWINAPI cublasShutdown (void);
}


using namespace Rcpp;

RcppExport SEXP gpuQR_magma(SEXP X_)
{    
    // Input
    NumericMatrix X(X_);

    // Initialize magma and cublas
    magma_init();
    cublasInit();
    
    // Declare variables 
    int info, lwork, n_rows = X.nrow(), n_cols = X.ncol(), min_mn = std::min(n_rows, n_cols);
    double tmp[1];
    NumericVector scale(min_mn);
    double *work; 
    // Query workspace size
    magma_int_t  nb;




    nb     = magma_get_dgeqrf_nb( n_rows, n_cols );
    Rcpp::Rcout << "in here " << " " << nb  << " .... " << std::endl;


    magma_dmalloc_cpu( &work,    min_mn );

    magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]), &(work[0]), -1, &info); 
    Rcpp::Rcout << "in here " <<  " .... " << std::endl;
    lwork = work[0];
    Rcpp::Rcout << "in here " << " " << lwork << " .... " << std::endl;
//    NumericVector work(lwork);

    // Run QR decomposition
    magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]), &(work[0]), lwork, &info);

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
    cublasShutdown();

    // Output  
    return wrap(X);
}



