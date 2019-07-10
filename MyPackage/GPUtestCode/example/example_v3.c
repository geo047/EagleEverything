// This is a simple standalone example. See README.txt
#include <iostream> 
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>




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
    magma_queue_destroy( queues2[0] );  \ 
    magma_queue_destroy( queues2[1] );  \
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


// ------------------------------------------------------------
int main( int argc, char** argv )
{


   // Initialize Magma and print GPU environment
   magma_init (); 

   if (DEBUG)
       magma_print_environment();

   magma_int_t device = 0;  // default value of device number
   magma_int_t ngpu = 4;

   // Create queues
   magma_int_t ndevices;
   magma_device_t devices[ ngpu ];
   magma_getdevices( devices, ngpu, &ndevices );

   if (ngpu > ndevices){
     std::cout << " Error: Number of GPUs is greater than the number of devices that have been detected. " << std::endl;
     return 0;
   }

   // save in environment variable, so magma_num_gpus() picks it up
   char env_num_gpus[20];  // space for "MAGMA_NUM_GPUS=", 4 digits, and nil
   setenv( "MAGMA_NUM_GPUS", env_num_gpus, true );

   #ifdef HAVE_CUBLAS
      magma_setdevice( device );
   #endif

   // queue for default device
    magma_queue_t   queue;
    magma_queue_t   queues2[3];  // 2 queues + 1 extra NULL entry to catch errors

   #ifdef HAVE_CUBLAS
      // handle for directly calling cublas
       cublasHandle_t  handle;
   #endif



   // create queues on this device
   // 2 queues + 1 extra NULL entry to catch errors
    magma_queue_create( devices[ device ], &queues2[ 0 ] );
    magma_queue_create( devices[ device ], &queues2[ 1 ] );
    queues2[ 2 ] = NULL;

    queue = queues2[ 0 ];

    #ifdef HAVE_CUBLAS
       // handle for directly calling cublas
       handle = magma_queue_get_cublas_handle( queue );
    #endif


// Variables 
double gpu_time , cpu_time ;
magma_int_t n=500 , n2=n*n;
magma_int_t  ldda = magma_roundup( n, 32 );

// Convert X into a form that Magma can deal with of type double *
// Have to do it in stages
//double const *xx = X.begin();  // doesnt work in muiltiple GPU world because contents have to get copied into xx
double  xx[n2];

// copy contents to CPU that can be changed 
for (int j=0; j< 500 ; j++){
  for (int i=0; i < 500 ; i++){
     xx[i+n*j] = rand()*100;
  }
}

double *r; //  r - nxn matrices on the host
double *h_r ;   // nxn matrix with eigenvector on host
double *d_r ; // nxn matrix on the device
double * h_work ; // workspace
magma_int_t lwork ; // h_work size
magma_int_t * iwork ; // workspace
magma_int_t liwork ; // iwork size
double *w1  ; // w1  - vectors of eigenvalues
double error , work [1]; // used in difference computations
magma_int_t ione = 1, info ;
double mione = -1.0;
magma_int_t incr = 1;
magma_int_t ISEED [4] = {0 ,0 ,0 ,1}; // seed
CHECK_MALLOC(magma_dmalloc_cpu (&w1 ,n)); // host memory for real
CHECK_MALLOC(magma_dmalloc_cpu (&r,n2 )); // host memory for r
CHECK_MALLOC(magma_dmalloc_cpu (&h_r,n2 )); // host memory for r
CHECK_MALLOC(magma_dmalloc_cpu (& d_r ,n2 )); // device memory for d_r
// ---> CHECK_MALLOC( magma_malloc( (void**) &d_r_array,    ndevices * sizeof(double*) ));

// Query for workspace sizes
double   aux_work [1];
magma_int_t aux_iwork [1];
   PRINT("in here");
CHECK_GPU(magma_dsyevd_m (ngpu,  
                          MagmaVec , 
                          MagmaLower ,
                          n, xx   ,
                          n,
                          w1,
                         aux_work ,
                         -1, 
                         aux_iwork,
                         -1,
                         &info))
   PRINT("in here");
   PRINT(info);
//CHECK_GPU(magma_dsyevd_m (ngpu,  MagmaVec , MagmaLower ,n,NULL ,
//                         n,NULL ,
//                         aux_work ,
//                         -1, aux_iwork ,-1,&info ));

magma_int_t nb = magma_get_dsytrd_nb(n) ;
PRINT(nb);





if (info != MAGMA_SUCCESS){
     // clean up and shut down 
     SHUTDOWN;
     return 0;
};

// aux_work[0]  =  3200240001;
lwork = ( magma_int_t ) aux_work [0];
PRINT( ( long long int ) aux_work[0]);
PRINT("------");

liwork = aux_iwork [0];
PRINT("liwork");
PRINT(liwork);
PRINT("------");
magma_imalloc_cpu( &iwork, liwork ); 
magma_dmalloc_cpu (& h_work , lwork ); // memory for workspace

PRINT("lwork");
PRINT(lwork);
PRINT("------");





// Doing eigenvalue calculation in GPU land

CHECK_GPU(magma_dsyevd_m (ngpu,  
                          MagmaVec , 
                          MagmaLower ,
                          n, xx   ,
                          n,
                          w1,
                         h_work ,
                         lwork, 
                         iwork,
                         liwork,
                         &info))





PRINT("after gpu calc");
PRINT(xx[0]);
PRINT(xx[1]);
PRINT(xx[2]);
PRINT(xx[3]);
PRINT(xx[4]);

//CHECK_GPU( magma_dsyevd_m(ngpu, MagmaVec,MagmaLower,n,d_r,
//                          n,w1,r,n,h_work,
// lwork,iwork,liwork,&info));
if (info != MAGMA_SUCCESS){
    // clean up and shut down 
     SHUTDOWN;
     return 0;
};
if (DEBUG){
     std::cout << " First 10 eigenvalues (descending order)... " << std::endl;
     for(int i=0; i<10; i++){
        std::cout << w1[i] << " " ;
     }
     std::cout << "\n\n" << std::endl;

     std::cout << " First 10 rows and columns of the eigenvector matrix (residing in GPU memory) " << std::endl;
  // --->    magma_dprint_gpu(10,10,d_r_array,n, queue);
};


// clean up and shut down 
SHUTDOWN;




// Return list structure 
return 0;

} 






