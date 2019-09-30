#include <Rcpp.h>
#include<magma_v2.h>
#include "magma_lapack.h"
#include<iostream>
#include<fstream>
#include<chrono>



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


  magma_int_t n, n2,  info = 0;  // define MAGMA_ILP64 to get these as 64 bit integers
  magma_vec_t jobv ;  // tell the function if we want eigenvectors or not
  if (wantvectors)
                jobv = MagmaVec;  // we want vectors returned, not just eigenvalues
        else
                jobv = MagmaNoVec ;

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

  double * h_values; 
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_values, n) )
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

auto start = std::chrono::high_resolution_clock::now();
size = bfile.tellg();
    std::cout << "size=" << size << "\n"; 
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
std::cout << "size = bfile.tellg() =  " << elapsed.count() << "\n";



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


 double *A;
 magma_dmalloc_cpu(&A, n2);   // to house data that will be lost 
 for(int i=0; i<n2; i++){
     A[i] = h_vectors[i];
 }



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



  // AWG
  // get work size
  magma_int_t lda, LDVL, LDVR;
            LDVL=n;
            LDVR=n;
            lda   = n;
            n2    = lda*n;
            double *wr;
            magma_dmalloc_cpu(&wr, n);
            double *wl;
            magma_dmalloc_cpu(&wl, n);
            double *VL;
            magma_dmalloc_cpu(&VL, LDVL*n);
            double *VR;
            magma_dmalloc_cpu(&VR, LDVR*n);
            magma_int_t lwork = -1;
            double work[1];


  // get optimal workspace size  
  magma_dgeev(MagmaNoVec, MagmaVec, n, NULL, lda, NULL, NULL, NULL, LDVL, NULL, LDVR, work, -1, &info );
std::cout << "Your kidding ... " << std::endl;
  
lwork = work[0];
double *h_work;
magma_dmalloc_cpu(&h_work, lwork);


// Perform analysis 
magma_dgeev(MagmaNoVec, MagmaVec, n, A , lda, h_values , wl, VL, LDVL, h_vectors, LDVR, h_work,
              lwork, &info);

std::cout << "Your kidding ... " << std::endl;










  return info ;


}



