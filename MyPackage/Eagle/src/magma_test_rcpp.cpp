#include<Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include "magma_lapack.h"
#include "magma_v2.h"



using namespace Rcpp;


//--------------------------------------------
// Magma test code 
//--------------------------------------------
// [[Rcpp::export]]
Eigen::MatrixXd magma_test_rcpp ( )
{
 magma_init();

} // end function 




