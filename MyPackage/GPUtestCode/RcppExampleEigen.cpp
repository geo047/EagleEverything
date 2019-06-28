#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>



using namespace Rcpp;

// Potential issue
// First, you are accessing the matrix row-wise, whereas MAGMA and LAPACK use column-wise ordering. That is,
//
// A[ i + j*size ]


// [[Rcpp::export]]
Rcpp::NumericVector   gpuEigen_magma(Rcpp::NumericMatrix X)
{


magma_init (); // initialize Magma
magma_queue_t queue = NULL ;
magma_int_t dev =0;
magma_queue_create (dev ,& queue );
double gpu_time , cpu_time ;
// magma_int_t n=8192 , n2=n*n;
magma_int_t n=X.nrow() , n2=n*n;
magma_int_t  ldda = magma_roundup( n, 32 );


double *a, *r; // a, r - nxn matrices on the host
double *d_r ; // nxn matrix on the device
double * h_work ; // workspace
magma_int_t lwork ; // h_work size
magma_int_t * iwork ; // workspace
magma_int_t liwork ; // iwork size
double *w1 , *w2; // w1 ,w2 - vectors of eigenvalues
double error , work [1]; // used in difference computations
magma_int_t ione = 1, info ;
double mione = -1.0;
magma_int_t incr = 1;
magma_int_t ISEED [4] = {0 ,0 ,0 ,1}; // seed
magma_dmalloc_cpu (&w1 ,n); // host memory for real
magma_dmalloc_cpu (&w2 ,n); // eigenvalues
magma_dmalloc_cpu (&a,n2 ); // host memory for a
magma_dmalloc_cpu (&r,n2 ); // host memory for r
// ---> magma_dmalloc (& d_r ,n2 ); // device memory for d_r
magma_dmalloc (& d_r ,ldda * n ); // device memory for d_r
// Query for workspace sizes
double aux_work [1];
magma_int_t aux_iwork [1];
// --->magma_dsyevd_gpu ( MagmaVec , MagmaLower ,n,d_r ,n,w1 ,r,n, aux_work ,
// --->.-1, aux_iwork ,-1,& info );
magma_dsyevd_gpu ( MagmaVec , MagmaLower ,n,d_r ,ldda ,w1 ,r,n, aux_work ,
-1, aux_iwork ,-1,& info );


std::cout << aux_work[0] << std::endl;
std::cout << aux_iwork[0] << std::endl;

lwork = ( magma_int_t ) aux_work [0];
liwork = aux_iwork [0];
iwork =( magma_int_t *) malloc ( liwork * sizeof ( magma_int_t ));
magma_dmalloc_cpu (& h_work , lwork ); // memory for workspace


// Create the matrix a and copy a -> r

// define a, r // [1 0 0 0 0 ...]
for(int i=0;i<n;i ++){ // [0 2 0 0 0 ...]
a[i*n+i]=( double )(i +1); // a = [0 0 3 0 0 ...]
r[i*n+i]=( double )(i +1); // [0 0 0 4 0 ...]
}



// Testing .....
magma_dsetmatrix ( n, n, r, n, d_r ,n, queue ); // copy a -> d_r
magma_dprint_gpu(5,5, d_r, n , queue);


// copy X to host memory 
double *h_TEST;
std::cout << "Testing copy idea  allocating h_TEST with dmalloc_cpu " << std::endl;
magma_dmalloc_cpu (&h_TEST,n2 ); // host memory for a
std::cout << "Testing copy idea " << std::endl;

//for(int i=0; i<n; i++){
//  for(int j=0; j<n; j++){
//     h_TEST[i+j*n] <- X[i,j];
//     std::cout << X[i,j] << std::endl;
//  }
//}
 magma_int_t lda = n;
 magma_int_t ldq = n;
std::cout << "Copying X to h_TEST " << std::endl;
 lapackf77_dlacpy( "Full", &n,   &n, &X[0], &lda, &h_TEST[0], &ldq );
std::cout << "Print h_TEST " << std::endl;
// // --->  magma_dprint(n,n,h_TEST,n);

std::cout << "Print X " << std::endl;
double const* XX = X.begin();
// --->  magma_dprint(n,n,XX ,n);



// -----> magma_dsetmatrix ( n, n, h_TEST , n, d_r ,n, queue ); // copy a -> d_r



// compute the eigenvalues and eigenvectors for a symmetric ,
// real nxn matrix ; Magma version
//magma_dsyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,
//lwork,iwork,liwork,&info);





//magma_dsetmatrix ( n, n, h_X, n, d_r ,n, queue ); // copy a -> d_r
 magma_dsyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,
 lwork,iwork,liwork,&info);
std::cout << "Wahoo !!!! " <<std::endl;


NumericVector ans =  NumericVector(w1,w1 + n  );
// ans.attr("dim") = Dimension(M, N);



magma_free_cpu (w1 ); // free host memory
magma_free_cpu (w2 ); // free host memory
magma_free_cpu (a); // free host memory
magma_free_cpu (r); // free host memory
magma_free_cpu ( h_work ); // free host memory
magma_free_cpu (h_TEST); // free host memory
magma_free (d_r ); // free device memory
magma_queue_destroy ( queue ); // destroy queue
magma_finalize (); // finalize Magma



return ans ;





} 
