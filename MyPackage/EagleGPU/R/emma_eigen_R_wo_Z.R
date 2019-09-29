
emma.eigen.R.wo.Z <-  function (K, X, ngpu=1)
{
    n <- nrow(X)
    q <- ncol(X)
    # S <- -1 *  X %*% solve(crossprod(X, X)) %*% t(X)
    S <- -1 *  X %*% magmaSolve(Xmat = (crossprod(X, X)), ngpu=ngpu, printInfo=FALSE ) %*% t(X)

    diag(S) <- diag(S) + 1
    diag(K) <- diag(K) + 1  ## old code: K + dn
    gc()
    if(ngpu > 0){
       #eig <- rcppMagmaSYEVD::eigen_mgpu(S %*% K  %*% S, symmetric = TRUE, only_values=FALSE)
       XX <-  (S %*% K  %*% S)

       eig <- magmaEigen(Xmat = XX   , ngpu=ngpu, printInfo=TRUE)
    } else {
       eig <- eigen(S %*% K %*% S, symmetric = TRUE)
    }


    stopifnot(!is.complex(eig$values))
    return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[,
        1:(n - q)]))
}


