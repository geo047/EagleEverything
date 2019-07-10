// This is a simple standalone example. See README.txt
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>


// set to 1 if debug information is desired
#define DEBUG 1



/***************************************************************************//**
 * Macros 
 */

#define PRINT(x)  \
     std::cout << x << std::endl



#ifdef HAVE_CUBLAS
#define SHUTDOWN   \
     magma_free_cpu (w1 ); \
     magma_free_cpu (r); \
     magma_free_cpu ( h_work ); \
     magma_free_cpu ( h_r ); \
     magma_free_cpu (d_r); \
     magma_finalize (); \
     magma_free_cpu(iwork) ; \
    queue = NULL;    \
    magma_queue_destroy( queues2[0] );  \ 
    magma_queue_destroy( queues2[1] );  \
    queues2[0] = NULL;     \
    queues2[1] = NULL;   \
    handle = NULL  
#else
   # define SHUTDOWN      \
     magma_free_cpu (w1 ); \
     magma_free_cpu (r); \
     magma_free_cpu ( h_work ); \
     magma_free_cpu ( h_r ); \
     magma_free_cpu (d_r); \
     magma_finalize (); \
     magma_free_cpu(iwork); \
    queue = NULL;    \
    magma_queue_destroy( queues2[0] );  \ 
    magma_queue_destroy( queues2[1] );  \
    queues2[0] = NULL;     \
    queues2[1] = NULL   
#endif






int main( int argc, char** argv )
{
  // Initialize Magma and print GPU environment
   magma_init ();

   if (DEBUG)
       magma_print_environment();

  magma_int_t ncpu = 16;
  magma_int_t ngpu = 1;

  magma_int_t *iwork;
  magma_int_t lda, info;
  magma_int_t n=500   , n2=n*n;
  magma_int_t  ldda = magma_roundup( n, 32 );
  magma_int_t threads = ncpu;
  magma_int_t  Nfound, lwork, liwork;
  double *h_A, *h_R, *h_work;
  double *w1, *w2;
  lda = n;
  magma_int_t il, iu;
  double vl, vu;

  // Convert X into a form that Magma can deal with of type double *
  // Have to do it in stages
  //double const *xx = X.begin();  // doesnt work in muiltiple GPU world because contents have to get copied into xx
  double  xx[n2];
  magma_dmalloc_cpu( &h_A,   n2 );

 /* Initialize the matrix */
//  magma_generate_matrix( opts, n, n, h_A, n );
 

  // copy contents to CPU that can be changed 
  for (int j=0; j< n; j++){
    for (int i=j; i < n; i++){
        h_A[i+n*j] = std::rand() ;
        h_A[j+n*i] = h_A[i+n*j];
    }
  }
  // changing diags 
  for (int j=0; j< n; j++){
        h_A[j+n*j] = std::rand();
    }




  // Get sizes 
  magma_dsyevdx_getworksize(n, threads, MagmaVec,
                                     &lwork,
                                     &liwork);






   /* Allocate host memory for the matrix */
   magma_dmalloc_cpu( &w1,    n );
   magma_dmalloc_cpu( &w2,    n );
   magma_imalloc_cpu( &iwork, liwork );

   magma_dmalloc_pinned( &h_R,    n2    );
   magma_dmalloc_pinned( &h_work, lwork );





   lapackf77_dlacpy( MagmaFullStr, &n, &n, h_A, &n, h_R, &n );


  if (ngpu == 1) {
                    //printf("calling dsyevdx_2stage 1 GPU\n");
                    magma_dsyevdx_2stage( MagmaVec, MagmaRangeAll, MagmaLower , n,
                                          h_R, lda,
                                          NULL, vu, il, NULL,
                                          &Nfound, w1,
                                          h_work, lwork,
                                          iwork, liwork,
                                          &info );
                } else {
                  PRINT("in here");
                    magma_dsyevdx_2stage_m( ngpu, MagmaNoVec, MagmaRangeAll, MagmaUpper , n,  
                                            h_R, lda,
                                            NULL , vu, il, iu,
                                            &Nfound, w1,
                                            h_work, lwork,
                                            iwork, liwork,
                                            &info );
                }

PRINT(info);
PRINT(Nfound);
PRINT(w1[0]);
PRINT(w1[1]);
PRINT(w1[2]);
PRINT(w1[3]);
PRINT("----------------");
  for (int j=0; j< n; j++){
    for (int i=0; i < n; i++){
//        PRINT(h_R[i+n*j]);
    }
  }
PRINT("----------------");





    return 0;



}
