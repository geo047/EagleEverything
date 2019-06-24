#include <Rcpp.h>
#include<magma_v2.h>
#include<magma_lapack.h>
#include<magma_operators.h>


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif



using namespace Rcpp;




// [[Rcpp::export]]
//RcppExport SEXP gpuQR_magma(SEXP X)
//Rcpp::NumericMatrix   gpuQR_magma(Rcpp::NumericMatrix X)
RcppExport SEXP   gpuQR_magma(SEXP X_)

{
 NumericMatrix X(X_);
   Rcpp::Rcout << X[0,0] << std::endl;
   Rcpp::Rcout << X[0,1] << std::endl;
   Rcpp::Rcout << X[0,2] << std::endl;
   Rcpp::Rcout << X[0,3] << std::endl;

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
// print matrix contents on host device
std::cout << "THIS HAS WORKED - Generated matrix via lapack thingie " << std::endl;
magma_dprint(M,N, h_A, lda);






 /* ====================================================================
               Performs operation using MAGMA
    =================================================================== */


// MAGMA
 double *d_a ; // d_a mxn matrix on the device
magma_dmalloc (&d_a , ldda*M); // device memory for d_a
magma_dsetmatrix ( M, N, h_A ,M,d_a ,ldda , queue ); // copy a -> d_a
// see if matrix on gpu matches matrix on host
std::cout << "THIS HAS WORKED - see if matrix on gpu matches matrix on host " << std::endl;
magma_dprint_gpu(M,N, d_a, ldda, queue);






// QR algorithm
magma_dgeqrf2_gpu( M, N, d_a, ldda, tau, &info );
std::cout << "Results coming out of dgeqrf2 " << std::endl;
 magma_dprint_gpu(M, N, d_a , ldda, queue);   // print contents of matrix on gpu device

 double *r;
  magma_dmalloc_pinned (&r,n2 ); // host memory for r
 magma_dgetmatrix ( M, N, d_a ,ldda ,r, M , queue ); // copy d_a -> r
std::cout << "printing out contents of r " << std::endl;
magma_dprint(M,N, r, lda);
  


std::cout << "About to start transfer or r into R land ... " << std::endl;
std::cout << n2 << std::endl;
std::cout << sizeof(r) << std::endl;
std::cout << sizeof(*r) << std::endl;

NumericVector ans =  NumericVector(r,r + n2  );
ans.attr("dim") = Dimension(M, N);

return ans;



}


