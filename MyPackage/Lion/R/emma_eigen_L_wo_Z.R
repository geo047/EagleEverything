emma.eigen.L.wo.Z <- function (K, ngpu=1)
{

   doMagmaEigen <- FALSE
   if (ngpu == 1) {
     if ( nrow(K) > 4000 )
           doMagmaEigen <- TRUE
   } else if (ngpu == 2) {
     if (nrow(K) > 4500 )
           doMagmaEigen <- TRUE
   } else if (ngpu == 3) {
     if (nrow(K) > 4500 )
           doMagmaEigen <- TRUE
  } else if (ngpu > 3) {
     if (nrow(K) > 6000)
           doMagmaEigen <- TRUE
  } else {
    message(" Error: the number of gpu needs to be set to a sensible number. \n")
    return(NULL)
 }

   if (doMagmaEigen){
    eig <- magmaEigen(Xmat=K, ngpu=ngpu, printInfo=TRUE)
   } else {
    eig <- eigen(K, symmetric = TRUE)
   }



    return(list(values = eig$values, vectors = eig$vectors))
}


