// This is a simple standalone example. See README.txt
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>










// ------------------------------------------------------------
int main( int argc, char** argv )
{

magma_init (); // initialize Magma
magma_queue_t queue = NULL ;
magma_int_t dev =0;
magma_queue_create (dev ,& queue );
double gpu_time , cpu_time ;
magma_int_t n=150 , n2=n*n;
float *a, *r; // a, r - nxn matrices on the host
float *d_r; // nxn matrix on the device
float * h_work ; // workspace
magma_int_t lwork ; // h_work size
magma_int_t * iwork ; // workspace
magma_int_t liwork ; // iwork size
float *w1 , *w2; // w1 ,w2 - vectors of eigenvalues
float error , work [1]; // used in difference computations
magma_int_t ione = 1, info ;
float mione = -1.0f;
magma_int_t incr = 1;
magma_int_t ISEED [4] = {0 ,0 ,0 ,1}; // seed
magma_smalloc_cpu (&w1 ,n); // host memory for real
magma_smalloc_cpu (&w2 ,n); // eigenvalues
magma_smalloc_cpu (&a,n2 ); // host memory for a
magma_smalloc_cpu (&r,n2 ); // host memory for r
magma_smalloc (& d_r ,n2 ); // device memory for d_r
// Query for workspace sizes
float aux_work [1];
magma_int_t aux_iwork [1];
magma_ssyevd_gpu ( MagmaVec , MagmaLower ,n,d_r ,n,w1 ,r,n, aux_work ,
-1, aux_iwork ,-1,& info );
lwork = ( magma_int_t ) aux_work [0];
liwork = aux_iwork [0];
iwork =( magma_int_t *) malloc ( liwork * sizeof ( magma_int_t ));
magma_smalloc_cpu (& h_work , lwork ); // memory for workspace
// Randomize the matrix a and copy a -> r
lapackf77_slarnv (& ione ,ISEED ,&n2 ,a);
lapackf77_slacpy ( MagmaFullStr ,&n ,&n,a ,&n,r ,&n);
magma_ssetmatrix ( n, n, a, n, d_r ,n, queue ); // copy a -> d_r
std::cout << "d_r " << std::endl;
magma_sprint_gpu(5,5,d_r,n,queue);

// compute the eigenvalues and eigenvectors for a symmetric ,
// real nxn matrix ; Magma version
gpu_time = magma_sync_wtime ( NULL );
magma_ssyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,
lwork,iwork,liwork,&info);

free (w1 ); // free host memory
free (w2 ); // free host memory
free (a); // free host memory
free (r); // free host memory
free ( h_work ); // free host memory
magma_free (d_r ); // free device memory
magma_queue_destroy ( queue ); // destroy queue
magma_finalize (); // finalize Magma



}