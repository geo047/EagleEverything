#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>



using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List   gpuEigen_magma(const Rcpp::NumericMatrix&  X)
{

// Initialize Magma and print GPU environment
magma_init (); 
magma_print_environment();


// Initialize the queue
magma_queue_t queue = NULL ;
magma_int_t dev =0;
magma_queue_create (dev ,& queue );


// Variables 
double gpu_time , cpu_time ;
magma_int_t n=X.nrow() , n2=n*n;
magma_int_t  ldda = magma_roundup( n, 32 );


double *r; //  r - nxn matrices on the host
double *h_r ;   // nxn matrix with eigenvector on host
double *d_r ; // nxn matrix on the device
double * h_work ; // workspace
magma_int_t lwork ; // h_work size
magma_int_t * iwork ; // workspace
magma_int_t liwork ; // iwork size
double *w1  ; // w1  - vectors of eigenvalues
double error , work [1]; // used in difference computations
magma_int_t ione = 1, info ;
double mione = -1.0;
magma_int_t incr = 1;
magma_int_t ISEED [4] = {0 ,0 ,0 ,1}; // seed
magma_dmalloc_cpu (&w1 ,n); // host memory for real
magma_dmalloc_cpu (&r,n2 ); // host memory for r
magma_dmalloc_cpu (&h_r,n2 ); // host memory for r
magma_dmalloc (& d_r ,n2 ); // device memory for d_r
//----> magma_dmalloc (& d_r ,ldda * n ); // device memory for d_r


// Query for workspace sizes
double aux_work [1];
magma_int_t aux_iwork [1];
magma_dsyevd_gpu ( MagmaVec , MagmaLower ,n,d_r ,n,w1 ,r,n, aux_work ,
   -1, aux_iwork ,-1,& info );
// -----> magma_dsyevd_gpu ( MagmaVec , MagmaLower ,n,d_r ,ldda ,w1 ,r,n, aux_work ,
// -----> -1, aux_iwork ,-1,& info );
lwork = ( magma_int_t ) aux_work [0];
liwork = aux_iwork [0];
iwork =( magma_int_t *) malloc ( liwork * sizeof ( magma_int_t ));
magma_dmalloc_cpu (& h_work , lwork ); // memory for workspace


// Convert X into a form that Magma can deal with of type double *
// Have to do it in stages
// double const *xx = &X(0,0);
double const *xx = X.begin();

// Copy contents of X that sit in CPU memory to GPU memory : xx -> d_r  
magma_dsetmatrix ( n, n, xx  , n, d_r ,n , queue ); 
std::cout << "Contents of d_r " << std::endl;
magma_dprint_gpu(5,5, d_r, n , queue);

// Doing eigenvalue calculation in GPU land
 magma_dsyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,
 lwork,iwork,liwork,&info);

std::cout << "EigenVectors on the GPU  " << std::endl;
magma_dprint_gpu(5,5,d_r,n, queue);




magma_dgetmatrix(n,n,d_r, n, h_r, n, queue);
std::cout << "EigenVectors on the CPU  " << std::endl;
magma_dprint(5,5, h_r,n );
std::cout << " -----------------------------  " << std::endl;



NumericVector vectors = NumericVector(h_r, h_r + n2);
std::cout << " in here " << std::endl;
vectors.attr("dim") = Dimension(n,n);

std::cout << " in here " << std::endl;
// Return list structure 
return Rcpp::List::create( 
       Rcpp::Named("values") =   NumericVector(w1,w1 + n  ) , 
       Rcpp::Named("vectors") =  vectors );

std::cout << " in here " << std::endl;

magma_free_cpu (w1 ); // free host memory
magma_free_cpu (r); // free host memory
magma_free_cpu ( h_work ); // free host memory
magma_free_cpu ( h_r ); // free host memory
magma_free (d_r ); // free device memory
magma_queue_destroy ( queue ); // destroy queue
magma_finalize (); // finalize Magma




} 

