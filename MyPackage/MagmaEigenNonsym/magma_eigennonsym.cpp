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


  // AWG
  // get work size
   magma_int_t n, n2;
  magma_int_t lda, LDVL, LDVR, info;
            n = numrows;
            LDVL=n;
            LDVR=n;
            lda   = n;
            n2    = lda*n;
            double *wr;
 if (  MAGMA_SUCCESS !=    magma_dmalloc_cpu(&wr, n) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for wr. Need more memory." << std::endl;
    return -1;
  }



            double *wl;
            if (  MAGMA_SUCCESS !=    magma_dmalloc_cpu(&wl, n) )
             {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for wl. Need more memory." << std::endl;
    return -1;
  }





//            double *VL;
//  if (  MAGMA_SUCCESS !=    magma_dmalloc_cpu(&VL, LDVL*n)  )
//  {
//    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for VL. Need more memory." << std::endl;
//    return -1;
//  }

            double *VR;
  if (  MAGMA_SUCCESS !=    magma_dmalloc_cpu(&VR, LDVR*n)  )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for VR. Need more memory." << std::endl;
    return -1;
  }

            magma_int_t lwork = -1;
            double work[1];






  magma_vec_t jobv ;  // tell the function if we want eigenvectors or not
  if (wantvectors)
                jobv = MagmaVec;  // we want vectors returned, not just eigenvalues
        else
                jobv = MagmaNoVec ;


  // get optimal workspace size  
  magma_dgeev(MagmaVec, MagmaNoVec, n, NULL, lda, NULL, NULL, NULL, LDVL, NULL, LDVR, work, -1, &info );       
           
  lwork = work[0];
  double *h_work;
 if (  MAGMA_SUCCESS !=  magma_dmalloc_cpu(&h_work, lwork))
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_work. Need more memory." << std::endl;
    return -1;
  }






 // ask for memory for eigenvalues
   double * h_values;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_values, n) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_values. Need more memory." << std::endl;
    return -1;
  }


  double * h_vectors;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_vectors, n2) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_vectors. Need more memory" << std::endl;
    return -1;
  }

  double * A;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&A, n2) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for A. Need more memory" << std::endl;
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

 size = n * n * sizeof(double);
  memblock = new char [size];
  bfile.seekg (0, std::ios::beg);
  bfile.read (memblock, size);
  bfile.close();

  A  =  (double*)memblock;  //reinterpret as doubles

  if (printInfo){
     std::cout << " First 5 rows and columns of the  X matrix " << std::endl;
     magma_dprint(5,5, A, n);
  }

  if(printInfo){
     std::cout << "Reading of X  matrix data complete.  " << std::endl;
  }





  if (numgpus == 1) {
     if (printInfo){
       std::cout << " About to begin single-gpu dsyevdx_2stage eigendecomposition calculation. " << std::endl;
     }
         magma_dgeev(MagmaNoVec, jobv , n, A, lda, h_values, wl, NULL, LDVL, h_vectors, LDVR, h_work,
              lwork, &info);

  } else {
    if (printInfo){
       std::cout << " About to begin multi-gpu dsyevdx_2stage eigendecomposition calculation. " << std::endl;
       std::cout << " Computation will use " << numgpus << " gpus. " << std::endl;
     }

    magma_dgeev_m( MagmaNoVec, jobv, n, A, lda, h_values, wl, NULL, LDVL , h_vectors, LDVR , h_work,
            lwork, &info ); 

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
  magma_free_cpu(h_work) ;
  magma_free_cpu(wr) ;
  magma_free_cpu(wl) ;
//  magma_free_cpu(VL) ;
  magma_free_cpu(VR) ;
  magma_free_cpu(A) ;

 magma_finalize();
  
 return 0;


}
                       
  



