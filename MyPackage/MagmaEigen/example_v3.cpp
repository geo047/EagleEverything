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


#ifdef HAVE_CUBLAS
#define SHUTDOWN   \
     magma_free_cpu (w1 ); \
     magma_free_cpu (r); \
     magma_free_cpu ( h_work ); \
     magma_free_cpu ( h_r ); \
     magma_free_cpu (d_r); \
     magma_finalize (); \
     magma_free_cpu(iwork) ; \
    queue = NULL;    \
    magma_queue_destroy( queues2[0] ); \  
    magma_queue_destroy( queues2[1] ); \ 
    queues2[0] = NULL;     \
    queues2[1] = NULL;   \
    handle = NULL  
#else
   # define SHUTDOWN      \
     magma_free_cpu (w1 ); \
     magma_free_cpu (r); \
     magma_free_cpu ( h_work ); \
     magma_free_cpu ( h_r ); \
     magma_free_cpu (d_r); \
     magma_finalize (); \
     magma_free_cpu(iwork); \
    queue = NULL;    \
    magma_queue_destroy( queues2[0] );  \ 
    magma_queue_destroy( queues2[1] );  \
    queues2[0] = NULL;     \
    queues2[1] = NULL   
#endif

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

 std::string X = argv[1];  // name of binary file with data
 int numrows = atoi(argv[2]); // number of rows 
 int numgpus = atoi(argv[3]); // number of rows 

 bool printInfo = (bool) atoi(argv[4]); // printInfo value (0 no print, 1 print)
 std::string fnamevec = argv[5];   
 std::string fnameval = argv[6];   
 bool wantvectors = (bool) atoi(argv[7]); 




 // Initialize Magma and print GPU environment
 magma_init (); 

  if (printInfo)
       magma_print_environment();

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
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_vectors. Need more memory" << std::endl;
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



// std::cout << n2 << std::endl;

 std::cout << h_vectors[250000-4] << std::endl;
 std::cout << h_vectors[250000-3] << std::endl;
 std::cout << h_vectors[250000-2] << std::endl;
 std::cout << h_vectors[250000-1] << std::endl;
 std::cout << h_vectors[250000] << std::endl;


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





  if ( MAGMA_SUCCESS !=  magma_dmalloc_cpu(&work, lwork))
  {
      std::cout   << " Error: magma_eigen() magma_malloc_cpu failed for work. Need more memory." << std::endl;
      return -1;
  }

  if ( MAGMA_SUCCESS !=  magma_imalloc_cpu(&iwork, liwork)){
      std::cout << " Error: magma_eigen() magma_malloc_cpu failed for iwork. Need more memory." << std::endl;
      return -1;
  }

        m = 0;

   double * h_values;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_values, n) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_values. Need more memory." << std::endl;
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



  magma_free_cpu (h_vectors)  ;
  magma_free_cpu (h_values)  ;
  magma_free_cpu(iwork)  ;
  magma_free_cpu(work) ;

  magma_queue_destroy ( queue );


  return info ;


}
                       
  



