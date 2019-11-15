
emma.eigen.R.wo.Z <-  function (K, X )
{
    n <- nrow(X)
    q <- ncol(X)
    # S <- -1 *  X %*% solve(crossprod(X, X)) %*% t(X)
     S <- -1 *  X %*%  chol2inv(chol(crossprod(X, X)))  %*% t(X)

    diag(S) <- diag(S) + 1
    diag(K) <- diag(K) + 1  ## old code: K + dn
    gc()

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
       XX <-  (S %*% K  %*% S)
       eig <- magmaEigen(Xmat = XX   ,  printInfo=FALSE)
   } else {
       eig <- eigen(S %*% K %*% S, symmetric = TRUE)
   }


    stopifnot(!is.complex(eig$values))
    return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[,
        1:(n - q)]))
}


