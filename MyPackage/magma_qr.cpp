#include <Rcpp.h>
#include<magma_v2.h>

#if defined(_OPENMP)
#include <omp.h>
#endif


/* Author: Andrew W. George
   Date:   16 Sep, 2019
   Purpose: to implement a multi-GPU version of the R function qr() which is very slow. 
*/




//--------------------------------------------
// Magma test code 
//--------------------------------------------
// [[Rcpp::export]]
 int  magma_qr(Rcpp::NumericMatrix  X , int numgpus, bool printInfo, std::string fname,  Rcpp::Function message )
{

  /*
  - Single- and mult-gpu  Magma code for performing QR factorisiation. 
  - A square data matrix is assumed. 
  - The data is read in from R and a copy is made otherwise the contents are overridden. 
  - The interface between R and Magma is 32 bits, even when Magma and R are built as 64 bit code. 
    Josh got around this by designing a server to run on shared memory. However, shared memory is limited to 
    matrices of around 46000 rows/cols. 
    I got around the problem by writing Q to disk, as a binary file, and reading this binary file back into R (very fast).   
    A simple solution to a very complicated problem. 
  */

   magma_init();


  if (printInfo){
   magma_print_environment();
  }

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
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_X, n2) )
  {
    message(" Error: magma_qr() magma_dmalloc_cpu failed for h_X. Need more memory.\n" );
    return(-1);
  }


 // Assign data to CPU 
  #if defined(_OPENMP)
     #pragma omp parallel for  
  #endif
  for (int j=0; j< n; j++){
    for (int i=0; i < n; i++){
        h_X[i+n*j] = X(i,j);
    }
  }




 if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu( &tau,     n  ) )
 {
    message(" Error: magma_qr() magma_dmalloc_cpu failed for tau. Need more memory. \n" );
    return(-1);

 }

 if(numgpus==1){

    if (  MAGMA_SUCCESS !=  magma_dmalloc( &dA,     ldda*n )  )
     {
          message(" Error: magma_qr() magma_malloc failed for dA. Most likely need more memory on GPU. Run on machine with multiple GPU. \n" );
          return -1;
     }


     if (  MAGMA_SUCCESS !=  magma_dmalloc( &dT,     ( 2*min_mn + magma_roundup( n, 32 ) )*nb )  )
     {
          message(" Error: magma_qr() magma_malloc failed for dT. Most likely need more memory on GPU. Run on machine with multiple GPU. \n" );
          return(-1);
     }


     magma_dsetmatrix(  n, n, h_X  , lda, dA, ldda, queue );
     if (printInfo) {
         message(" Single GPU \n ");
         message(" First 5 rows/columns of R data that is now sitting on the GPU \n ");
         magma_dprint_gpu(5,5, dA, n , queue);
     }


     magma_dgeqrf_gpu( n, n, dA, ldda, tau, dT, &info );

     if(info < 0){
        message(" magma_dgeqrf_gpu has returned error message ",   info, "\n");
        message(" The number identifies which parameter is the cause of the error. \n");
        return info;
     }



     if (printInfo ) {
         message(" Output from running magma_dgeqrf_gpu: contents still on GPU \n ");
         magma_dprint_gpu(5,5, dA, n , queue);
     }

     magma_dorgqr_gpu( n, n, n, dA, ldda, tau, dT, nb, &info );

     if(info < 0){
        message( " magma_dorgqr_gpu has returned error message \n" );
        message(" The number identifies which parameter is the cause of the error. \n");
        return info;
     }

     if (printInfo ) {
         message(" Q matrix. From running magma_dorgqr_gpu: contents still on GPU \n ");
         magma_dprint_gpu(5,5, dA, n , queue);
     }


     magma_dgetmatrix(  n, n ,dA,    ldda,  h_X  , lda, queue );
     if (printInfo ) {
         message( " Q matrix. GPU -> CPU. This is what is going to be passed into R land (via binary file). \n" );
         magma_dprint(5,5, h_X , n );
     }

     // Write Q matrix to disc in binary format
     if (printInfo){
         message("About to start writing to the tempoarary binary file. \n ");
     }
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

    if(info < 0){
        message( " magma_dgeqrf has returned error message " , info, "\n" );
        message( " The number identifies which parameter is the cause of the error. \n" );
        return info;
     }


     lwork = tmp[0];


    double * h_work ;

    if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&h_work, lwork) )
     {
        message(" Error: magma_qr() magma_dmalloc_cpu failed for h_work. Need more memory. \n" );
        return(0);
     }



   magma_dgeqrf_m(numgpus  , n, n, h_X   , n   , tau, h_work, lwork, &info );
   if(info < 0){
        message( " magma_dgeqrf_m has returned error message " ,info, "\n" );
        message( " The number identifies which parameter is the cause of the error. \n " );
        return info;
    }

   if (printInfo ) {
         message( " Output from running magma_dgeqrf_m: contents on CPU \n "  );
         magma_dprint(5, 5, h_X , n);

   }

   double *cmat;
    if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&cmat, n2) )
     {
        message(" Error: magma_qr() magma_dmalloc_cpu failed for cmat. Need more memory. \n" );
        return(0);
     }


   // Form identify matrix - used for forming Q matrix in next step
  #if defined(_OPENMP)
     #pragma omp parallel for  
  #endif
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

   if(info < 0){
        message( " magma_dormqr_m  has returned error message " , info, "\n" );
        message( " The number identifies which parameter is the cause of the error. \n" );
        return info;
    }


   if (printInfo ) {
         message( " Q matrix  from running magma_dormqr_m: contents on CPU \n"  );
          magma_dprint(5,5, cmat , n);
   }


   // Write Q matrix to disc in binary format
     if (printInfo){
         message( "About to start writing to the tempoarary binary file. \n " );
     }

   FILE* file = fopen(fname.c_str() , "wb");
    fwrite(&cmat[0], 1 , n*n*sizeof(double) , file);
    fclose(file);


   magma_free_cpu  (h_X);
   magma_free_cpu  (h_work);
   magma_free_cpu  (tau);
   magma_free_cpu  (cmat);




}


  magma_queue_destroy ( queue );






  return (info) ;


}























