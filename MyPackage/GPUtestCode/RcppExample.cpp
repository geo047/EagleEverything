#include <Rcpp.h>
#include<magma_v2.h>
#include<magma_lapack.h>
#include<magma_operators.h>

using namespace Rcpp;




// [[Rcpp::export]]
Rcpp::NumericVector   gpuQR_magma(Rcpp::NumericMatrix X)
{

   // convert to C++ type
   double const* d_X = X.begin();

    // Initialize magma and cublas
    magma_init();

   magma_print_environment();

   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );
    const double c_zero    = MAGMA_D_ZERO;
   
    double  *tau,  tmp[1], unused[1];
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

    magma_dmalloc( &d_A,    ldda*N );

    size = (2*M + magma_roundup( N, 32 ) )*nb;
    magma_dmalloc( &dT, size );
    magmablas_dlaset( MagmaFull, size, 1, c_zero, c_zero, dT, size, queue );
 
 /* ====================================================================
               Performs operation using MAGMA
    =================================================================== */


// MAGMA
 double *d_a ; // d_a mxn matrix on the device
magma_dmalloc (&d_a , ldda*M); // device memory for d_a

magma_dsetmatrix ( M, N, d_X ,M,d_a ,ldda , queue ); // copy X -> d_a
std::cout << " This is what is now sitting in d_a on the GPU ready for analysis" << std::endl;
magma_dprint_gpu(M, N, d_a , ldda, queue);   // print contents of matrix on gpu device



// QR algorithm
magma_dgeqrf2_gpu( M, N, d_a, ldda, tau, &info );
std::cout << "Results coming out of dgeqrf2 - has d_a changed? " << std::endl;
magma_dprint_gpu(M, N, d_a , ldda, queue);   // print contents of matrix on gpu device

double *r;
magma_dmalloc_pinned (&r,n2 ); // host memory for r
magma_dgetmatrix ( M, N, d_a ,ldda ,r, M , queue ); // copy d_a -> r
std::cout << "Contents of r after moving contents d_a to r which is in host memory " << std::endl;
magma_dprint(M,N, r, lda);
  



NumericVector ans =  NumericVector(r,r + n2  );
ans.attr("dim") = Dimension(M, N);

return ans;



}


