// This is a simple standalone example. See README.txt

#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>
#include <cstdlib>
#include<fstream>
#include<chrono>





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

#define CHECK_GPU( err )   \
   if (err != MAGMA_SUCCESS ){                         \
       std::cout <<"\n" << std::endl;    \
       std::cout << "\nError:    " << std::endl;    \
       std::cout << "  Eigenvalue calculation has failed. " << std::endl;  \
       std::cout << "  Things to try to solve this problem are:            " << std::endl;  \
       std::cout << "    1. Reduce the number of individuals to see if this is a GPU memory issue. " << std::endl; \
       std::cout << "    2. Ensure that you do not have perfectly correlated fixed effects in the model. " << std::endl; \
       std::cout << "       This will cause collinearity between columns in the model matrix. \n\n " << std::endl;     \
   }



/*
 Command line parameters
--------------------------
std::string  X  
long  numrows, 
int numgpus, 
bool printInfo, 
std::string fnamevec,  
std::string fnameval,  
bool wantvectors
*/



int main( int argc, char** argv )
{
auto start = std::chrono::high_resolution_clock::now();

 // Input received from command line
 std::string X = argv[1];  // name of binary file with data
 int numrows = atoi(argv[2]); // number of rows 
 int numgpus = atoi(argv[3]); // number of gpu 

 bool printInfo = (bool) atoi(argv[4]); // printInfo value (0 no print, 1 print)
 std::string fnameval = argv[5];  // name of binary file containing eigenvector answers 
 std::string fnamevec = argv[6];  // name of binrary file containing eigenvalues.  
 bool wantvectors = (bool) atoi(argv[7]);    // flag - if eigenvectors are wanted.  



 // Initialize Magma and print GPU environment
 magma_init (); 

 if (printInfo)
       magma_print_environment();



  magma_int_t n, n2, lwork, info = 0;  // define MAGMA_ILP64 to get these as 64 bit integers
  magma_vec_t jobv ;  // tell the function if we want eigenvectors or not
  if (wantvectors)
                jobv = MagmaVec;  // we want vectors returned, not just eigenvalues
        else
                jobv = MagmaNoVec ;

  double *work ;
  magma_int_t liwork, *iwork, m;
  double vl = 0.0, vu = 0.0, abstol = 0.0;


   n = numrows ;   // this is X.rows() but passed in from command line
   n2     = n*n;
   magma_int_t il = 1;
   magma_int_t iu = n;

  double * h_vectors;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_vectors, n2) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_vectors. Need more memory" << std::endl;
    return -1;
  }

  // X matrix being passed in as a binary file
  std::streampos size;
  char * memblock;

  std::ifstream bfile ( X.c_str() , std::ios::in|std::ios::binary|std::ios::ate);


//size = bfile.tellg();
//    std::cout << "size=" << size << "\n";









  if(printInfo){
     std::cout << "About to begin reading in the X matrix data ... " << std::endl;
  }

//auto start = std::chrono::high_resolution_clock::now();
  size = n * n * sizeof(double);
  memblock = new char [size];
  bfile.seekg (0, std::ios::beg);
  bfile.read (memblock, size);
  bfile.close();

  h_vectors =  (double*)memblock;  //reinterpret as doubles

//auto finish = std::chrono::high_resolution_clock::now();
//std::chrono::duration<double> elapsed = finish - start;
//std::cout << "Read X matrix binary file =   " << elapsed.count() << "\n";

  if(printInfo){
     std::cout << "Reading of X  matrix data complete.  " << std::endl;
  }

  if (printInfo){
     std::cout << " First 5 rows and columns of the  X matrix " << std::endl;
     magma_dprint(5,5, h_vectors, n);
  }



  // ask for optimal size of work arrays 
  lwork = -1; liwork = -1;
  double aux_work[1];
  magma_int_t aux_iwork[1];

//start = std::chrono::high_resolution_clock::now();
  magma_dsyevdx_2stage( jobv, MagmaRangeAll, MagmaLower , n, NULL  , n,
                      vl, vu, il, iu, &m, NULL , aux_work, -1, aux_iwork, -1 , &info);

//finish = std::chrono::high_resolution_clock::now();
//elapsed = finish - start;
//std::cout << "Workspace query =    " << elapsed.count() << "\n";



   lwork = ( magma_int_t ) aux_work [0];
   liwork = aux_iwork [0];


  // ask for cpu memory for work variables
  if ( MAGMA_SUCCESS !=  magma_dmalloc_cpu(&work, lwork))
  {
      std::cout   << " Error: magma_eigen() magma_malloc_cpu failed for work. Need more memory." << std::endl;
      return -1;
  }

  if ( MAGMA_SUCCESS !=  magma_imalloc_cpu(&iwork, liwork)){
      std::cout << " Error: magma_eigen() magma_malloc_cpu failed for iwork. Need more memory." << std::endl;
      return -1;
  }

  // ask for memory for eigenvalues
   double * h_values;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_values, n) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_values. Need more memory." << std::endl;
    return -1;
  }




  if (numgpus == 1) {
     if (printInfo){
       std::cout << " About to begin single-gpu dsyevdx_2stage eigendecomposition calculation. " << std::endl;
     }
//start = std::chrono::high_resolution_clock::now();
      magma_dsyevdx_2stage( jobv, MagmaRangeAll, MagmaLower , n,
                     h_vectors   , n, vl, vu, il, iu, &m, h_values, work, lwork,
                      iwork, liwork, &info);
//finish = std::chrono::high_resolution_clock::now();
//elapsed = finish - start;
//std::cout << "1 GPU: magma_dsyevdx_2stage =    " << elapsed.count() << "\n";



  } else {
     if (printInfo){
       std::cout << " About to begin multi-gpu dsyevdx_2stage eigendecomposition calculation. " << std::endl;
       std::cout << " Computation will use " << numgpus << " gpus. " << std::endl;
     }

//start = std::chrono::high_resolution_clock::now();
      magma_dsyevdx_2stage_m(numgpus, jobv, MagmaRangeAll , MagmaLower , n,
                      h_vectors, n, vl, vu, il, iu, &m, h_values,
                      work, lwork, iwork, liwork, &info);
//finish = std::chrono::high_resolution_clock::now();
//elapsed = finish - start;
//std::cout << "m GPU: magma_dsyevdx_2stage_m =    " << elapsed.count() << "\n";



  }

  if (printInfo){
      std::cout << "The largest 4 eigenvalues found"<< std::endl ;
      for (int t1 = 1 ; t1 < 5 ; t1++)
      {
          std::cout << h_values[n-t1] << ",\t" ;
       }
       std::cout << std::endl ;
   }


   if( printInfo){
       std::cout << " About to write eigenvalues to binary file. " << std::endl;
    }


  if (fnameval.empty()){
    std::cout << " Error: in magma_eigen. Name of binary file for the eigenvalues has not been specified." << std::endl;
    return -1;
  }


    FILE* file = fopen(fnameval.c_str() , "wb");
    if (file == NULL){
       std::cout << " Error: in magma_eigen. Binary file for eigenvalues has failed to open. " << std::endl;
       return -1;
    } else {
       fwrite(&h_values[0], 1 , n*sizeof(double) , file);
       fclose(file);
    }



   if( printInfo){
       std::cout << " About to write eigenvectors  to binary file. " << std::endl;
    }


  if(wantvectors){

    if (fnamevec.empty()){
      std::cout << " Error: in magma_eigen. Name of binary file for the eigenvectors has not been specified. \n" << std::endl;
      return -1;
    }

//start = std::chrono::high_resolution_clock::now();
    FILE* file = fopen(fnamevec.c_str() , "wb");
    if (file == NULL){
       std::cout << " Error: in magma_eigen. Binary file for eigenvectors has failed to open." << std::endl;
       return -1;
    } else {
       fwrite(&h_vectors[0], 1 , n*n*sizeof(double) , file);
       fclose(file);
    }
  }
 

  // tidying up 
  magma_free_cpu (h_vectors)  ;
  magma_free_cpu (h_values)  ;
  magma_free_cpu(iwork)  ;
  magma_free_cpu(work) ;
//  delete [] aux_work;
//  delete [] aux_iwork;

 magma_finalize();

auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
std::cout << "Whole thing  =   " << elapsed.count() << "\n";
  
 return 0;


}
                       
  



