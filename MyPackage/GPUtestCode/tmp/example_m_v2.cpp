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

#define CHECK_MALLOC( err )                                                                        \
    if (err != MAGMA_SUCCESS){                                                                     \
        std::cout << "\nError: Eigenvalue calculation has failed due to failure of memory allocation.\n" << std::endl;  \
     }



#define SHUTDOWN   \
     magma_free_cpu ( work ); \
     magma_free_cpu ( h_R ); \
     magma_free_cpu (tau); \
     magma_queue_destroy ( queue ); \
     magma_finalize (); 





int main( int argc, char** argv )
{


   // convert to C++ type
  


  // Initialize magma and cublas
   magma_init();
   magma_print_environment();

   // Initialize the queue
   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );


   double  *work, *h_R, *h_A, *tau ;
   magma_int_t M, N, n2, lda, ldda,  info,  nb;
   magma_int_t ngpu=4;

    M = 1000;
    N = 1000;
   lda    = M;
   n2     = lda*N;
   ldda = ((M +31)/32)*32; // ldda = m if 32 divides m
   nb = magma_get_dsytrd_nb(N); 

   magma_int_t lwork =  std::max( 2*N + N*nb, 1 + 6*N + 2*N*N );
   // mallocs
   magma_dmalloc_cpu(&work, lwork);
   magma_dmalloc_cpu(&tau, N);
   magma_dmalloc_cpu(&h_R, n2);
  magma_dmalloc_cpu( &h_A,    N*lda  );

 // Assign data to CPU 
//  for (int j=0; j< N; j++){
//    for (int i=0; i < N; i++){
//         h_R[i+N*j] = X(i,j);
//    }
//  }

 // Assign data to CPU 
  for (int j=0; j< N; j++){
    for (int i=j; i < N; i++){
         h_R[i+N*j] = (double) (std::rand() % 1000) ;
         h_R[j+N*i] = h_R[i+N*j];
    }
  }

 magma_dprint(10,10, h_R, N); 




  /* ====================================================================
               Performs operation using MAGMA
  =================================================================== */


  magma_int_t *iwork;

  magma_int_t  Nfound,  liwork;
  double *h_work;
  double *w1, *w2;
  magma_int_t il=0, iu=0;
  double vl=0, vu=0;
  magma_dmalloc_cpu( &w1,    N );   // eigenvalues found [out]
   magma_dmalloc_cpu( &w2,    N );
  liwork  = 3 + 5*N; 
magma_imalloc_cpu( &iwork, liwork );
magma_dmalloc_cpu(&h_work, lwork);
// ---> magma_imalloc( &iwork, liwork );
// ---> magma_dmalloc(&h_work, lwork);

std::cout << "in here ... " << std::endl;


   magma_dsyevdx_m( ngpu, MagmaVec, MagmaRangeAll, MagmaLower , N,
                                            h_R, 
                                            lda,
                                            vl , vu, il, iu,
                                            &Nfound, w1,
                                            h_work, lwork,
                                            iwork, liwork,
                                            &info );

std::cout << "Wow ... " << std::endl;



  if (info != 0) {
       printf("magma_dsyevdx_m returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
       }

    

for(int i=0; i< 3; i++){
  for(int j=0; j< 3; j++){
    std::cout << h_R[i+N*j] << std::endl;
  }
}

// clean up and shut down 
 SHUTDOWN;




// Return list structure 
return 0;




}
