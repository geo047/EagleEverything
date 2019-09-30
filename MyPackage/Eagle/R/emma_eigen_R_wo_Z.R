
emma.eigen.R.wo.Z <-  function (K, X, ngpu=0)
{
    n <- nrow(X)
    q <- ncol(X)
    dn <- diag(n)
    start <- Sys.time()
    S <- dn - X %*% solve(crossprod(X, X)) %*% t(X)
    end <- Sys.time()
    cat(" Solve : ", end - start, "\n")
    gc()
#    if(ngpu > 0){
#     if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
#       eig <- rcppMagmaSYEVD::eigen_mgpu(S %*% (K + dn) %*% S, symmetric = TRUE, only_values=FALSE)
#     }
#    } else {
    start <- Sys.time()
       eig <- eigen(S %*% (K + dn) %*% S, symmetric = TRUE)
    end <- Sys.time()
    cat(" Eigen calc  : ", end - start, "\n")
#    }


    stopifnot(!is.complex(eig$values))
    return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[,
        1:(n - q)]))
}


