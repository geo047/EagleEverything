#include<Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include "magma_lapack.h"
#include "magma_v2.h"



using namespace Rcpp;


//--------------------------------------------
// Magma test code 
//--------------------------------------------
// [[Rcpp::export]]
int  magma_qr ( const Rcpp::NumericMatrix  X , int numgpus=1, bool printInfo=false, std::string fname )
{

    // convert to C++ typ, std::string fnamee
   double const* h_X = X.begin();   // this is a column-wise vector which magma likes

   magma_init();
   magma_print_environment();

   // Initialize the queue
   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );


  magma_int_t n, n2, lwork, info = 0;  // define MAGMA_ILP64 to get these as 64 bit integers
  double   *T  ;

  double  *tau;
  magmaDouble_ptr dA,  dT;

 n = X.nrow() ;

  magma_int_t lda, ldda, min_mn, nb;
  lda  = n;
  ldda = n;  // multiple of 32 by default
  n2 = lda*n;
  min_mn = n;
  nb = magma_get_dgeqrf_nb( n, n );
  lwork  = n*nb;

 std::cout << "In here " <<std::endl;
 if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu( &tau,     n  ) )
 {
    shrd_server->error_and_die(" MAGMA_QR_SERVER Error: magma_qr_mgpu() magma_malloc_cpu failed. Most likely need more memory." );
 }

 if(numgpus==1){

    if (  MAGMA_SUCCESS !=  magma_dmalloc( &dA,     ldda*n )  )
     {
          shrd_server->error_and_die(" MAGMA_QR_SERVER Error: magma_qr_mgpu() magma_malloc failed. Most likely need more memory on GPU. Run on machine with multiple GPU." );
     }

 std::cout << "In here " <<std::endl;

     if (  MAGMA_SUCCESS !=  magma_dmalloc( &dT,     ( 2*min_mn + magma_roundup( n, 32 ) )*nb )  )
     {
    shrd_server->error_and_die(" MAGMA_QR_SERVER Error: magma_qr_mgpu() magma_malloc failed. Most likely need more GPU memory. Run on machine with multiple GPU." );
     }

 std::cout << "about to magma_dsetmatrix  " <<std::endl;

     magma_dsetmatrix(  n, n, h_X  , lda, dA, ldda, queue );
     if (printInfo) {
         std::cout << " Single GPU "  << std::endl;
         std::cout << " First 5 rows/columns of R data that is now sitting on the GPU "  << std::endl;
         magma_dprint_gpu(5,5, dA, n , queue);
     }


     magma_dgeqrf_gpu( n, n, dA, ldda, tau, dT, &info );
     if (printInfo ) {
         std::cout << " Output from running magma_dgeqrf_gpu: contents still on GPU "  << std::endl;
         magma_dprint_gpu(5,5, dA, n , queue);
     }

     magma_dorgqr_gpu( n, n, n, dA, ldda, tau, dT, nb, &info );

     if (printInfo ) {
         std::cout << " Q matrix. From running magma_dorgqr_gpu: contents still on GPU "  << std::endl;
         magma_dprint_gpu(5,5, dA, n , queue);
     }


     magma_dgetmatrix(  n, n ,dA,    ldda,  h_X  , lda, queue );
     if (printInfo ) {
         std::cout << " Q matrix. GPU -> CPU. This is what is being passed into R land."  << std::endl;
         magma_dprint(5,5, h_X , n );
     }

    magma_free ( dT )  ;
    magma_free ( dA )  ;
    magma_free_cpu ( tau )  ;

 } else {
     // query for workspace size
     lwork = -1;
     double tmp[1];
     magma_dgeqrf( n, n, NULL, n, NULL, tmp, lwork, &info );
     lwork = tmp[0];


    double * h_work ;

    if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&h_work, lwork) )
     {
    shrd_server->error_and_die(" MAGMA_QR_SERVER Error: magma_qr_mgpu() magma_dmalloc_cpu failed. Most likely need more CPU memory." );
     }


   magma_dgeqrf_m(shrd_server->_numgpus  , n, n, h_X   , n   , tau, h_work, lwork, &info );
   if (printInfo ) {
         std::cout << " Output from running magma_dgeqrf_m: contents on CPU "  << std::endl;
         magma_dprint(5, 5, h_X , n);

   }

   double *cmat;
    if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&cmat, n2) )
     {
        shrd_server->error_and_die(" MAGMA_QR_SERVER Error: magma_qr_mgpu() magma_dmalloc_cpu failed. Most likely need more CPU memory." );
     }


   // Form identify matrix - used for forming Q matrix in next step
   for( magma_int_t col=0; col < n; col++){
      for( magma_int_t row=0; row< n; row++){
      if( row == col){
          cmat[row + n * col] = 1;
      } else {
          cmat[row + n * col] = 0;

      }
     }
   }


   magma_dormqr_m  ( shrd_server->_numgpus , MagmaLeft, MagmaNoTrans, n, n, n, h_X, n, tau, cmat , n, h_work, lwork, &info ) ;

   if (printInfo ) {
         std::cout << " Q matrix  from running magma_dormqr_m: contents on CPU "  << std::endl;
         magma_dprint(5,5, cmat , n);
   }

  // Copy Q matrix from cmat to h_X (R land)
lapackf77_dlacpy( MagmaFullStr ,&n ,&n,cmat  ,&n, h_X ,&n); // a- >r


   magma_free_cpu  (h_work);
   magma_free_cpu  (tau);






}


  magma_queue_destroy ( queue );

// Write Q matrix to disc in binary format

std::cout << "Starting writing fwrite ... " << std::endl;
FILE* file = fopen(fname.c_str() , "wb");
    fwrite(&h_R[0], 1 , N*N*sizeof(double) , file);
    fclose(file);
auto endTime2 = std::chrono::high_resolution_clock::now();





  return (info) ;


}



















} // end function 




