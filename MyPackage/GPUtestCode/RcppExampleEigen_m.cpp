#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>

// set to 1 if debug information is desired
#define DEBUG 1


using namespace Rcpp;

/***************************************************************************//**
 * Macros 
 */

#define PRINT(x)  \
     std::cout << x << std::endl

#define CHECK_MALLOC( err )                                                                        \
    if (err != MAGMA_SUCCESS){                                                                     \
        Rcpp::Rcout << "\nError: Eigenvalue calculation has failed due to failure of memory allocation.\n" << std::endl;  \
     }

#define CHECK_GPU( err )   \
   if (err != MAGMA_SUCCESS ){                         \
       Rcpp::Rcout <<"\n" << std::endl;    \
       Rcpp::Rcout << "\nError:    " << std::endl;    \
       Rcpp::Rcout << "  Eigenvalue calculation has failed. " << std::endl;  \
       Rcpp::Rcout << "  Things to try to solve this problem are:            " << std::endl;  \
       Rcpp::Rcout << "    1. Reduce the number of individuals to see if this is a GPU memory issue. " << std::endl; \
       Rcpp::Rcout << "    2. Ensure that you do not have perfectly correlated fixed effects in the model. " << std::endl; \
       Rcpp::Rcout << "       This will cause collinearity between columns in the model matrix. \n\n " << std::endl;     \
   }


#define SHUTDOWN   \
     magma_free_cpu (w1 ); \
     magma_free_cpu (w2 ); \
     magma_free_cpu (h_A); \
     magma_free_cpu ( iwork ); \
     magma_free_cpu ( h_R ); \
     magma_free_cpu (h_work); \
     magma_finalize (); 





// [[Rcpp::export]]
Rcpp::List   gpuEigen_magma(const Rcpp::NumericMatrix&  X)
{
  // Initialize Magma and print GPU environment
   magma_init ();

   if (DEBUG)
       magma_print_environment();

  magma_int_t ncpu = 16;
  magma_int_t ngpu = 4;

  magma_int_t *iwork;
  magma_int_t lda, info;
  magma_int_t n=X.nrow()  , n2=n*n;
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


//  for (int j=0; j< n; j++){
//    for (int i=0; i < n; i++){
//        h_A[i+n*j] =  X(i,j);
//    }
//  }



  // Get sizes 
  magma_dsyevdx_getworksize(n, threads, MagmaVec,
                                     &lwork,
                                     &liwork);






   /* Allocate host memory for the matrix */
   magma_dmalloc_cpu( &w1,    n );
   magma_dmalloc_cpu( &w2,    n );
   magma_imalloc_cpu( &iwork, liwork );

   // ---> magma_dmalloc_pinned( &h_R,    n2    );
   magma_dmalloc_cpu( &h_R,    n2    );
   magma_dmalloc_cpu( &h_work, lwork );





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
                    magma_dsyevdx_2stage_m( ngpu, MagmaVec, MagmaRangeAll, MagmaLower , n,  
                                            h_R, lda,
                                            NULL , vu, il, iu,
                                            &Nfound, w1,
                                            h_work, lwork,
                                            iwork, liwork,
                                            &info );
                }

PRINT(info);
PRINT(w1[0]);
PRINT(w1[1]);
PRINT(w1[2]);
PRINT(w1[3]);
//PRINT("----------------");
//  for (int j=0; j< n; j++){
//    for (int i=0; i < n; i++){
//        PRINT(h_R[i+n*j]);
//    }
//  }
PRINT("----------------");



SHUTDOWN;


    return 0;



}
