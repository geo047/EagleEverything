#include <Rcpp.h>
#include<magma_v2.h>
//#include<magma_lapack.h>
//#include<magma_operators.h>




/*

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "magma_v2.h"
// #include "magma_myown.h"


*/

 // int  magma_qr(Eigen::Map<Eigen::MatrixXd>  X , int numgpus, bool printInfo, std::string fname,  Rcpp::Function message )


//--------------------------------------------
// Magma test code 
//--------------------------------------------
// [[Rcpp::export]]
 int  magma_qr(Rcpp::NumericMatrix  X , int numgpus, bool printInfo, std::string fname,  Rcpp::Function message )
{
    // convert to C++ typ, std::string fnamee
  // ----> double * h_X = X.begin();
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

 n = X.rows() ;

  magma_int_t lda, ldda, min_mn, nb;
  lda  = n;
  ldda = n;  // multiple of 32 by default
  n2 = lda*n;
  min_mn = n;
  nb = magma_get_dgeqrf_nb( n, n );
  lwork  = n*nb;





  double * h_X;
  magma_dmalloc_cpu(&h_X, n2);


 // Assign data to CPU 
  for (int j=0; j< n; j++){
    for (int i=0; i < n; i++){
        h_X[i+n*j] = X(i,j);
    }
  }




 if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu( &tau,     n  ) )
 {
    message(" Error: magma_qr() magma_dmalloc_cpu failed. Need more memory." );
    return(0);

 }

 if(numgpus==1){

    if (  MAGMA_SUCCESS !=  magma_dmalloc( &dA,     ldda*n )  )
     {
          message(" Error: magma_qr() magma_malloc failed. Most likely need more memory on GPU. Run on machine with multiple GPU." );
          return(0);
     }

 std::cout << "In here " <<std::endl;

     if (  MAGMA_SUCCESS !=  magma_dmalloc( &dT,     ( 2*min_mn + magma_roundup( n, 32 ) )*nb )  )
     {
          message(" Error: magma_qr() magma_malloc failed. Most likely need more memory on GPU. Run on machine with multiple GPU." );
          return(0);
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

     // Write Q matrix to disc in binary format
     std::cout << "Starting writing fwrite ... " << std::endl;
     FILE* file = fopen(fname.c_str() , "wb");
     fwrite(&h_X[0], 1 , n*n*sizeof(double) , file);
     fclose(file);

    magma_free ( dT )  ;
    magma_free ( dA )  ;
    magma_free_cpu ( tau )  ;
    magma_free_cpu( h_X);
 } else {
     // query for workspace size
     lwork = -1;
     double tmp[1];
     magma_dgeqrf( n, n, NULL, n, NULL, tmp, lwork, &info );
     lwork = tmp[0];


    double * h_work ;

    if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&h_work, lwork) )
     {
        message(" Error: magma_qr() magma_dmalloc_cpu failed. Need more memory." );
        return(0);
     }

   magma_dprint(n,n, h_X, n);


   magma_dgeqrf_m(numgpus  , n, n, h_X   , n   , tau, h_work, lwork, &info );
   if (printInfo ) {
         std::cout << " Output from running magma_dgeqrf_m: contents on CPU "  << std::endl;
         magma_dprint(5, 5, h_X , n);

   }

   double *cmat;
    if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&cmat, n2) )
     {
        message(" Error: magma_qr() magma_dmalloc_cpu failed. Need more memory." );
        return(0);
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


   magma_dormqr_m  ( numgpus , MagmaLeft, MagmaNoTrans, n, n, n, h_X, n, tau, cmat , n, h_work, lwork, &info ) ;

   if (printInfo ) {
         std::cout << " Q matrix  from running magma_dormqr_m: contents on CPU "  << std::endl;
          magma_dprint(5,5, cmat , n);
   }


   // Write Q matrix to disc in binary format
   std::cout << "Starting writing fwrite ... " << std::endl;
   FILE* file = fopen(fname.c_str() , "wb");
    fwrite(&cmat[0], 1 , n*n*sizeof(double) , file);
    fclose(file);

   // write out tail end of cmat for checking - hard coded to be 60000x60000
//  for( magma_int_t col=59996; col < 60000; col++){
//      for( magma_int_t row=59996; row< 60000; row++){
//          std::cout << cmat[row + n * col] << " " ;
//     }
//      std::cout << std::endl;
//   }



   magma_free_cpu  (h_X);
   magma_free_cpu  (h_work);
   magma_free_cpu  (tau);
   magma_free_cpu  (cmat);




}


  magma_queue_destroy ( queue );






  return (info) ;


}























