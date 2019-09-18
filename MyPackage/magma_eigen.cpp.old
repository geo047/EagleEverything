#include <Rcpp.h>
#include<magma_v2.h>
#include "magma_lapack.h"
#include<iostream>
#include<fstream>

#if defined(_OPENMP)
#include <omp.h>
#endif

/* Author: Andrew W. George
   Date:   16 Sep, 2019
   Purpose: to implement a multi-GPU version of the R function qr() which is very slow. 
*/


//  int  magma_eigen(Rcpp::NumericMatrix  X , int numgpus, bool printInfo, std::string fnamevec,  std::string fnameval,  Rcpp::Function message , bool wantvectors)


//--------------------------------------------
// Magma test code 
//--------------------------------------------
// [[Rcpp::export]]
 int  magma_eigen(std::string  X , long  numrows, int numgpus, bool printInfo, std::string fnamevec,  std::string fnameval,  Rcpp::Function message , bool wantvectors)
{


  magma_init();


  if (printInfo){
   magma_print_environment();
  }

   // Initialize the queue
   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );


  magma_int_t n, n2, lwork, info = 0;  // define MAGMA_ILP64 to get these as 64 bit integers
  magma_vec_t jobv ;  // tell the function if we want eigenvectors or not
  if (wantvectors)
                jobv = MagmaVec;  // we want vectors returned, not just eigenvalues
        else
                jobv = MagmaNoVec ;

  double *work ;
  magma_int_t liwork, *iwork, m;
  double vl = 0.0, vu = 0.0, abstol = 0.0;


   // n = X.rows() ;
   n = numrows ;
   n2     = n*n;
   magma_int_t il = 1;
   magma_int_t iu = n;

  double * h_vectors; 
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_vectors, n2) )
  {
    message(" Error: magma_eigen() magma_dmalloc_cpu failed for h_vectors. Need more memory.\n" );
    return -1;
  }

/*

 // Assign data to CPU 
  for (int j=0; j< n; j++){
    for (int i=0; i < n; i++){
        h_vectors[i+n*j] = X(i,j);
    }
  }
*/

std::streampos size;
char * memblock;

std::ifstream bfile ( X.c_str() , std::ios::in|std::ios::binary|std::ios::ate);


size = bfile.tellg();
    std::cout << "size=" << size << "\n"; 

memblock = new char [size];
bfile.seekg (0, std::ios::beg);
bfile.read (memblock, size);
bfile.close();

std::cout << "the entire file content is in memory \n";
//double* double_values = (double*)memblock;//reinterpret as doubles
//
//for(int i=0; i<=10; i++)
//{
//double value = double_values[i];
//std::cout << "value ("<<i<<")=" << value << "\n";
//}
h_vectors =  (double*)memblock;  //reinterpret as doubles





// FILE *pBinFile;
// size_t result;
// pBinFile = fopen ( X.c_str()   , "rb" );
// result  = fread (h_vectors ,n2, 1    ,pBinFile);
// std::cout << n2 << std::endl;

 std::cout << h_vectors[250000-4] << std::endl;
 std::cout << h_vectors[250000-3] << std::endl;
 std::cout << h_vectors[250000-2] << std::endl;
 std::cout << h_vectors[250000-1] << std::endl;
 std::cout << h_vectors[250000] << std::endl;



// if (result != n2) {
//    message(" Error: magma_eigen() fread of matrix data has failed. \n");
//    return -1;
// }
// std::cout << "result = " << result << std::endl;
  std::cout << "This should e the contents of the X matrix " << std::endl;
  magma_dprint(5,5, h_vectors, n);


  magma_int_t threads = 1;
  #if defined(_OPENMP)
        #pragma omp parallel
        {
            threads = omp_get_num_threads();
        }
  #endif



  // ask for optimal size of work arrays 
  lwork = -1; liwork = -1;
//  magma_dsyevdx_getworksize(n, threads, wantvectors, &lwork, &liwork);

double aux_work[1];
magma_int_t aux_iwork[1];

magma_dsyevdx_2stage( jobv, MagmaRangeAll, MagmaLower , n,
                     NULL  , n,
                      vl, vu, il, iu,
                      &m, NULL ,
                      aux_work, -1,
                      aux_iwork, -1 ,
                      &info);

lwork = ( magma_int_t ) aux_work [0];
liwork = aux_iwork [0];





  if ( MAGMA_SUCCESS !=  magma_dmalloc(&work, lwork))
  {
      message(" Error: magma_eigen() magma_malloc_cpu failed for work. Need more memory.\n" );
      return -1;
  }

  if ( MAGMA_SUCCESS !=  magma_imalloc(&iwork, liwork)){
      message(" Error: magma_eigen() magma_malloc_cpu failed for iwork. Need more memory.\n" );
      return -1;
  }

        m = 0;

   double * h_values;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_values, n) )
  {
    message(" Error: magma_eigen() magma_dmalloc_cpu failed for h_values. Need more memory.\n" );
    return -1;
  }

 


  if (numgpus == 1) {
    std::cout << "IN HERE" << std::endl;



//double  *r;
//magma_dmalloc_cpu (&r,n2 );
//lapackf77_dlacpy ( MagmaFullStr ,&n ,&n, xx   ,&n,  r ,&n);
//magma_dprint(5,5, r, n);


      magma_dsyevdx_2stage( jobv, MagmaRangeAll, MagmaLower , n,
                     h_vectors   , n,
                      vl, vu, il, iu,
                      &m, h_values,
                      work, lwork,
                      iwork, liwork,
                      &info);
       std::cout << "What is going on ... " << info << std::endl;
  } else {
      //printf("calling dsyevdx_2stage_m %ld GPU\n", (long int) opts.ngpu);
      magma_dsyevdx_2stage_m(numgpus, jobv, MagmaRangeAll , MagmaUpper , n,
                      h_vectors, n,
                      vl, vu, il, iu,
                      &m, h_values,
                      work, lwork,
                      iwork, liwork,
                      &info);
  }

  std::cout << "The largest 4 eigenvalues found"<< std::endl ;
  for (int t1 = 1 ; t1 < 5 ; t1++)
  {
      std::cout << h_values[n-t1] << ",\t" ;
  }
  std::cout << std::endl ;


  if (fnameval.empty()){
    message(" Error: in magma_eigen. Name of binary file for the eigenvalues has not been specified. \n");
    return -1;
  }


    FILE* file = fopen(fnameval.c_str() , "wb");
    if (file == NULL){
       message(" Error: in magma_eigen. Binary file for eigenvalues has failed to open. \n");
       return -1;
    } else {
       fwrite(&h_values[0], 1 , n*sizeof(double) , file);
       fclose(file);
    }

  if(wantvectors){

    if (fnamevec.empty()){
      message(" Error: in magma_eigen. Name of binary file for the eigenvectors has not been specified. \n");
      return -1;
    }

    FILE* file = fopen(fnamevec.c_str() , "wb");
    if (file == NULL){
       message(" Error: in magma_eigen. Binary file for eigenvectors has failed to open. \n");
       return -1;
    } else {
       fwrite(&h_vectors[0], 1 , n*n*sizeof(double) , file);
       fclose(file);
    }
  } 



  magma_free_cpu (h_vectors)  ;
  magma_free_cpu (h_values)  ;
  magma_free_cpu(iwork)  ;
  magma_free_cpu(work) ;

  magma_queue_destroy ( queue );


  return info ;


}



