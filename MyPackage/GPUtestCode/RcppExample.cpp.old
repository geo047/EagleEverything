#include<Rcpp.h>
#include<magma_v2.h>
#include<magma_lapack.h>
#include<magma_operators.h>

using namespace Rcpp;

 // RcppExport SEXP gpuQR_magma(SEXP X_)
 RcppExport NumericMatrix  gpuQR_magma(NumericMatrix X)
{    
    // Input
//    NumericMatrix X(X_);

    // Initialize magma and cublas
    magma_init();
 //   cublasInit();

   magma_print_environment();

   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );
    const double             d_neg_one = MAGMA_D_NEG_ONE;
    const double             d_one     = MAGMA_D_ONE;
    const double c_neg_one = MAGMA_D_NEG_ONE;
    const double c_one     = MAGMA_D_ONE;
    const double c_zero    = MAGMA_D_ZERO;
    const magma_int_t        ione      = 1;
   
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    double           Anorm, error=0, error2=0;
    double *h_A, *h_R, *tau, *h_work, tmp[1], unused[1];
    magmaDouble_ptr d_A, dT;
    magma_int_t M, N, n2, lda, ldda, lwork, info, min_mn, nb, size;
   
    M = X.rows();
    N = X.cols();
    min_mn = std::min(M, N); 
    lda    = M;
    n2     = lda*N;
    ldda = ((M +31)/32)*32; // ldda = m if 32 divides m
    nb     = magma_get_dgeqrf_nb( M, N );

    // query for workspace size
    lwork = -1;
    lapackf77_dgeqrf( &M, &N, unused, &M, unused, tmp, &lwork, &info );
    lwork = (magma_int_t)MAGMA_D_REAL( tmp[0] );

    magma_dmalloc_cpu( &tau,    min_mn );
    magma_dmalloc_cpu( &h_A,    n2     );
    magma_dmalloc_cpu( &h_work, lwork  );

    magma_dmalloc_pinned( &h_R,    n2     );

    magma_dmalloc( &d_A,    ldda*N );

    size = (2*M + magma_roundup( N, 32 ) )*nb;
    magma_dmalloc( &dT, size );
    magmablas_dlaset( MagmaFull, size, 1, c_zero, c_zero, dT, size, queue );
 
   /* Initialize the matrix */
  // hard to get this to work ... 
//   magma_generate_matrix( opts, M, N, h_A, lda );

// Randomize the matrix  Alternate mehtod for generating a matrix
magma_int_t ISEED [4] = {0 ,0 ,0 ,1}; // seed
lapackf77_dlarnv ( &ione , ISEED , &n2 , h_A ); // randomize a

 /* ====================================================================
               Performs operation using MAGMA
    =================================================================== */
            nb = magma_get_dgeqrf_nb( M, N );


// MAGMA
 double *d_a ; // d_a mxn matrix on the device
 magma_dmalloc (&d_a , ldda*M); // device memory for d_a
magma_dsetmatrix ( M, N, h_A ,M,d_a ,ldda , queue ); // copy a -> d_a
 magma_dgeqrf2_gpu( M, N, d_a, ldda, tau, &info );
Rcpp::Rcout << tau[0] << std::endl;
Rcpp::Rcout << tau[1] << std::endl;
Rcpp::Rcout << tau[2] << std::endl;

double *r;
magma_dmalloc_pinned (&r,n2 ); // host memory for r
magma_dgetmatrix ( M, N, d_a ,ldda ,r, M, queue ); // copy d_a -> r

Rcpp::NumericMatrix ans = wrap(r);
    return ans  ;

//
//
//
//    // Declare variables 
//    int info, lwork, n_rows = X.nrow(), n_cols = X.ncol(), min_mn = n_rows;
//    double tmp[1];
//    double  work[1];
//    NumericVector scale(min_mn);
//
//   // Xtra
//   magma_queue_t queue = NULL ;
//   magma_int_t dev =0;
//  magma_int_t ione = 1, lhwork ; // lhwork - workspace size
//   magma_queue_create (dev ,& queue );
//   double *d_a ; // d_a mxn matrix on the device
//   magma_int_t ldda ;
//   ldda = ((n_rows +31)/32)*32; // ldda = m if 32 divides m
//   double *tau, unused[1]  ; // scalars defining the elementary reflectors
//   double *hwork ;  // hwork - workspace ; tmp -used in
//
//
//// query for workspace size
//lwork = -1;
//lapackf77_dgeqrf( &n_rows, &n_cols, unused, &n_rows, unused, tmp, &lwork, &info );
//lwork = (magma_int_t)MAGMA_D_REAL( tmp[0] );
//std::cout << "lapackf77_dgeqrf"  << std::endl;
//
//   magma_dmalloc_cpu (& hwork , lhwork ); // host memory for hwork
//   magma_dmalloc_cpu (&tau , min_mn ); // host memory for tau
//   magma_dmalloc (&d_a , ldda*n_rows); // device memory for d_a
//
//
//
//std::cout << "ldda " << ldda << std::endl;
//   magma_dsetmatrix ( n_rows, n_cols, as<double *>(X) , n_rows, d_a ,ldda , queue ); // copy a -> d_a
//std::cout << "ldda " << ldda << std::endl;
//   std::cout << d_a[0] << std::endl;
//   std::cout << d_a[1] << std::endl;
//   std::cout << d_a[2] << std::endl;
//   std::cout << d_a[1] << std::endl;
//
//    // Query workspace size
//
//
//
//
//    // magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]), &(work[0]), -1, &info); 
//    //  magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(tau[0]), &(work[0]), -1, &info); 
//    magma_dgeqrf2_gpu( n_rows, n_cols, d_a , ldda, tau , &info);
//    Rcpp::Rcout << info << std::endl;
//
//    lwork = work[0];
//    NumericVector workN(lwork);
//
//
//
//
//    // Run QR decomposition
//    //             m       n      A        lda     tau           work       lwork   info
////     magma_dgeqrf(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]), &(workN[0]), lwork, &info);
////     magma_dmalloc_cpu (& scale , min_mn ); // host memory for tau
////  magma_dgeqrf2_gpu(n_rows, n_cols, &(X[0]), n_rows, &(scale[0]),  &info);
//
////    // Scale factor result    
////    for(int ii = 1; ii < n_rows; ii++)
////    {
////        for(int jj = 0; jj < n_cols; jj++)
////        {
////            // if(ii > jj) { X[ii + jj * n_rows] *= scale[jj]; 
////            if(ii > jj) { X[ii + jj * n_rows] *= tau[jj]; }
////          
////        }
////    }
//
//    // Scale factor result    
//    for(int ii = 1; ii < n_rows; ii++)
//    {
//        for(int jj = 0; jj < n_cols; jj++)
//        {
//            // if(ii > jj) { X[ii + jj * n_rows] *= scale[jj]; 
//            if(ii > jj) { 
//                Rcpp::Rcout << d_a[ii + jj * n_rows]   << std::endl;
//               X[ii + jj * n_rows] *= tau[jj]; 
//                Rcpp::Rcout << d_a[ii + jj * n_rows]   << std::endl;
//                Rcpp::Rcout << tau[jj]    << std::endl;
//                std::cout << "--------------" << std::endl;
//            }
//          
//        }
//    }
//
//    // Shutdown magma and cublas    
//    magma_finalize();
//   //  cublasShutdown();
//
//     // Output  
//    return wrap(X);
}


