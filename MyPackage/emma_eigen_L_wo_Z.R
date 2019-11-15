emma.eigen.L.wo.Z <- function (K)
{

   doMagmaEigen <- FALSE
   if (computer$ngpu == 1) {
     if ( nrow(K) > 4000 )
           doMagmaEigen <- TRUE
   } else if (computer$ngpu == 2) {
     if (nrow(K) > 4500 )
           doMagmaEigen <- TRUE
   } else if (computer$ngpu == 3) {
     if (nrow(K) > 4500 )
           doMagmaEigen <- TRUE
  } else if (computer$ngpu > 3) {
     if (nrow(K) > 6000)
           doMagmaEigen <- TRUE
  } else {
    message(" Error: the number of gpu needs to be set to a sensible number. \n")
    return(NULL)
 }

   if (doMagmaEigen){
    eig <- magmaEigen(Xmat=K,  printInfo=TRUE)
   } else {
    eig <- eigen(K, symmetric = TRUE)
   }



    return(list(values = eig$values, vectors = eig$vectors))
}


