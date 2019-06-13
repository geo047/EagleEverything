// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "readblock.h"

// includes, project
#include "magma_v2.h"
#include "magma_lapack.h"






//--------------------------------------------
// Calculation of transformed blup a values
//--------------------------------------------
// [[Rcpp::export]]
Eigen::MatrixXd magma_test_rcpp ( )
{
 magma_init();
 magma_print_environment();

} // end function 




