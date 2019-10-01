
emma.eigen.R.w.Z <-  function (Z, K, X, complete = TRUE, ngpu=0 )
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    n <- nrow(Z)
    t <- ncol(Z)
    q <- ncol(X)
    # matrix too large for GPU based code
     #SZ <- Z - X %*%   solve(crossprod(X, X))    %*% crossprod(X, Z)
     SZ <- Z - X %*%    chol2inv(chol(crossprod(X, X)))    %*% crossprod(X, Z)

    #eig <- eigen(K %*% crossprod(Z, SZ), symmetric = FALSE, EISPACK = TRUE)
    eig <- magmaEigenNonsym(Xmat = (K %*% crossprod(Z, SZ) ) , ngpu=ngpu, printInfo=FALSE)
    if (is.complex(eig$values)) {
        eig$values <- Re(eig$values)
        eig$vectors <- Re(eig$vectors)
    }
    #qr.X <- qr.Q(qr(X))
    qr.X <- magmaQR(Xmat= X, ngpu=ngpu, printInfo=FALSE)
    values = eig$values[1:(t - q)]
    #vectors = qr.Q(qr(cbind(SZ %*% eig$vectors[, 1:(t - q)], qr.X)), complete = TRUE)[,
    #    c(1:(t - q), (t + 1):n)]))
    vectors = magmaQR(Xmat= (cbind(SZ %*% eig$vectors[, 1:(t - q)], qr.X)), ngpu=ngpu, printInfo=FALSE)[, c(1:(t - q), (t + 1):n)]
    return(list(values=values, vectors=vectors))

}


