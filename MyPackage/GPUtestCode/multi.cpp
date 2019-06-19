#include "Rcpp.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_operators.h"
#include "testings.h"


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
   // Using testing_dsyevd_gpu.cpp as a template for this code


    // Input
    NumericMatrix X(X_);
    double abstol; 
             
/* Local variables */
    double  *h_A, *h_R, *h_Z, *h_work, aux_work[1], unused[1];
    magma_int_t *iwork, *isuppz, *ifail, aux_iwork[1];
    double *w1, *w2, result[4]={0, 0, 0, 0}, eps,  runused[1];
     magmaDouble_ptr d_R, d_Z;

    // Initialize magma and cublas
    magma_init();
    magma_print_environment();    


    // Declare variables 
    int N = X.nrow(), n_cols = X.ncol(), min_mn ;
    magma_int_t  Nfound, info, lwork, liwork, lda, ldda;
    min_mn = N;   // because symmetric matrix

    NumericVector scale(min_mn);
    double *work; 
    // Query workspace size
    magma_int_t  nb;

Rcpp::Rcout << "in here ... " << std::endl;
      magma_opts_t test;
      Rcpp::Rcout << test << std::endl; 

     magma_opts opts;
     // setting argc and **argv manually here
     // based on 
     int numargs=5;
     char* args[] = { "prg_name", "--version", "1", "-L", "-JN" };
     opts.parse_opts( numargs, args );




     lda = N;
     ldda = magma_roundup( N, opts.align );  // multiple of 32 by default
     abstol = 0;  // auto, in zheevr
Rcpp::Rcout << "in here ... " << std::endl;

    magma_range_t range;
    magma_int_t il, iu;
    double vl, vu;
Rcpp::Rcout << "in here ... " << std::endl;
    
    opts.get_range( N, &range, &vl, &vu, &il, &iu );
Rcpp::Rcout << "in here ... " << std::endl;


   // query for workspace sizes

   magma_dsyevd_gpu( opts.jobz, opts.uplo,
                                  N, NULL, ldda, NULL,  // A, w
                                  NULL, lda,            // host A
                                  aux_work,  -1,
                                  aux_iwork, -1,
                                  &info );
   lwork  = (magma_int_t) MAGMA_D_REAL( aux_work[0] );
   liwork = aux_iwork[0];


   /* Allocate host memory for the matrix */
    magma_dmalloc_cpu( &h_A,    N*lda  );
    magma_dmalloc_cpu( &w1,     N      );
    magma_dmalloc_cpu( &w2,     N      );
    magma_imalloc_cpu( &iwork,  liwork );

    magma_dmalloc_pinned( &h_R,    N*lda  );
    magma_dmalloc_pinned( &h_work, lwork  );
    magma_dmalloc( &d_R,    N*ldda );

   /* Initialize the matrix */
    Rcpp::Rcout << "almost " <<  " .... " << std::endl;
    // matrix data is in h_A
   magma_generate_matrix( opts, N, N, h_A, lda );
   magma_dsetmatrix( N, N, h_A, lda, d_R, ldda, opts.queue );
    Rcpp::Rcout << "almost " <<  " .... " << std::endl;

 real_Double_t   gpu_time ; 
  gpu_time = magma_wtime();

  std::cout << opts.uplo << std::endl;
  magma_dsyevd_gpu( opts.jobz, opts.uplo,
                                  N, d_R, ldda, w1,
                                  h_R, lda,
                                  h_work, lwork,
                                  iwork, liwork,
                                  &info );


//    NumericVector work(lwork);


    // Scale factor result    
    for(int ii = 1; ii < N; ii++)
    {
        for(int jj = 0; jj < n_cols; jj++)
        {
            if(ii > jj) { X[ii + jj * N] *= scale[jj]; }
        }
    }

    // Shutdown magma and cublas    
    magma_finalize();
    cublasShutdown();

    // Output  
   //  return wrap(X);
}



