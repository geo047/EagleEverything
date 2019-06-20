#include "Rcpp.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include <cuda.h> 
#include <cuda_runtime.h> 

#if defined(__cplusplus)
extern "C" {
#endif

/* CUBLAS data types */
#define cublasStatus cublasStatus_t

cublasStatus CUBLASWINAPI cublasInit (void);
cublasStatus CUBLASWINAPI cublasShutdown (void);
}






using namespace Rcpp;

RcppExport SEXP gpuQR_magma(SEXP X_)
{    
   // Using testing_dsyevd_gpu.cpp as a template for this code


    // Input
    NumericMatrix X(X_);
    double abstol; 
           
     magma_init ();    // initialize Magma
     magma_queue_t queue=NULL; 
     magma_int_t dev=0; 
     magma_queue_create(dev,&queue);
// magma_int_t n=8192, n2=n*n;
magma_int_t n=1024, n2=n*n;
 float *a, *r;   // a, r - nxn matrices in magma_ssyevd_gpu
 float *a1;  // nxn matrix, copy of a used 
// double *a1;  // nxn matrix, copy of a used 
// float *h_work;   // workspace
double *h_work;   // workspace
magma_int_t lwork;   // h_work size
magma_int_t *iwork;   // workspace 
magma_int_t liwork;    // iwork size 
float *w1, *w2;   // w1, w2 -  vectors of eigenvalues
float *d_r;    // n x n matrix on the device
//double *w1, *w2;   // w1, w2 -  vectors of eigenvalues
double error , work [1];  // used in difference computations
magma_int_t ione = 1, i, j, info; 
float mione = -1.0;

magma_int_t ISEED[4] = {0,0,0,1};
magma_smalloc_cpu(&w1,n);
magma_smalloc_cpu(&w2,n);
magma_smalloc_cpu(&a,n2);
magma_smalloc_cpu(&r,n2);
magma_smalloc(&d_r,n2);

float   aux_work [1];
magma_int_t aux_iwork [1]; 

magma_ssyevd_gpu(MagmaVec,
                 MagmaLower,
                  n,
                  d_r,
                  n,
                  w1,
                  r,
                  n,
                  aux_work,
                   -1,
                  aux_iwork ,
                  -1,
                 &info ); 








}  
