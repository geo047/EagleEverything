#include <Rcpp.h>
#include<magma_v2.h>
#include<magma_lapack.h>
#include<magma_operators.h>

using namespace Rcpp;

// Potential issue
// First, you are accessing the matrix row-wise, whereas MAGMA and LAPACK use column-wise ordering. That is,
//
// A[ i + j*size ]


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
     magma_free_cpu ( work ); \
     magma_free_cpu ( h_R ); \
     magma_free_cpu (tau); \
     magma_queue_destroy ( queue ); \
     magma_finalize (); 





// [[Rcpp::export]]
Rcpp::List   gpuEigen_magma(const Rcpp::NumericMatrix  X)
{


   // convert to C++ type
// --->   double const* d_X = X.begin();   // this is a column-wise vector which magma likes
  


  // Initialize magma and cublas
   magma_init();
   magma_print_environment();

   // Initialize the queue
   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );


   double  *work, *h_R, *tau ;
   magma_int_t M, N, n2, lda, ldda,  info,  nb;
   magma_int_t ngpu=4;

    M = X.rows();
    N = X.cols();
   lda    = M;
   n2     = lda*N;
   ldda = ((M +31)/32)*32; // ldda = m if 32 divides m
   nb = magma_get_dsytrd_nb(N); 

   magma_int_t lwork =  std::max( 2*N + N*nb, 1 + 6*N + 2*N*N );

   // mallocs
   magma_dmalloc_cpu(&work, lwork);
   magma_dmalloc_cpu(&tau, N);
   magma_dmalloc_cpu(&h_R, n2);


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



// lapackf77_dlacpy ( MagmaFullStr ,&N ,&N,h_R ,&N, h_yourkidding  ,&N);
 magma_dprint(10,10, h_R, N); 


  /* ====================================================================
               Performs operation using MAGMA
  =================================================================== */

//  magma_dgeqrf_m( ngpu, M, N, h_R , N   , tau, work, lwork, &info );

  magma_int_t *iwork;

  magma_int_t  Nfound,  liwork;
  double *h_work;
  double *w1, *w2;
  magma_int_t il=0, iu=0;
  double vl=0, vu=0;
  magma_dmalloc_cpu(&h_work, lwork);
  magma_dmalloc_cpu( &w1,    N );   // eigenvalues found [out]
   magma_dmalloc_cpu( &w2,    N );
  liwork  = 3 + 5*N; 
   magma_imalloc_cpu( &iwork, liwork );



   magma_dsyevdx_m( ngpu, MagmaVec, MagmaRangeAll, MagmaLower , N,
                                            h_R, 
                                            lda,
                                            vl , vu, il, iu,
                                            &Nfound, w1,
                                            h_work, lwork,
                                            iwork, liwork,
                                            &info );




  if (info != 0) {
       printf("magma_dsyevdx_m returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
       }

    
// create R list structure
NumericVector values =  NumericVector(tau,tau + N  );
NumericVector vectors = NumericVector(h_R, h_R + n2);
vectors.attr("dim") = Dimension(N,N);


for(int i=0; i< 3; i++){
  for(int j=0; j< 3; j++){
    std::cout << h_R[i+N*j] << std::endl;
  }
}

// clean up and shut down 
 SHUTDOWN;




// Return list structure 
return Rcpp::List::create(
       Rcpp::Named("tau") =   values,
      Rcpp::Named("A") =  vectors );




}

