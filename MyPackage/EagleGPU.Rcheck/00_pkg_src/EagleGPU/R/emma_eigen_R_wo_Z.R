
emma.eigen.R.wo.Z <-  function (K, X, ngpu=0)
{
    n <- nrow(X)
    q <- ncol(X)
    #dn <- diag(n)
    #S <- dn - X %*% solve(crossprod(X, X)) %*% t(X)
    # memory efficient approach 
    S <- -1 *  X %*% solve(crossprod(X, X)) %*% t(X)
    diag(S) <- diag(S) + 1
    diag(K) <- diag(K) + 1  ## old code: K + dn
    gc()
    if(ngpu > 0){
     if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
       #eig <- rcppMagmaSYEVD::eigen_mgpu(S %*% (K + dn) %*% S, symmetric = TRUE, only_values=FALSE)
       eig <- rcppMagmaSYEVD::eigen_mgpu(S %*% K  %*% S, symmetric = TRUE, only_values=FALSE)
     }
    } else {
       #eig <- eigen(S %*% (K + dn) %*% S, symmetric = TRUE)
       eig <- eigen(S %*% K %*% S, symmetric = TRUE)
    }


    stopifnot(!is.complex(eig$values))
    return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[,
        1:(n - q)]))
}


