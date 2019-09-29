 emma.eigen.L.w.Z <- function (Z, K, complete = TRUE, ngpu=1)
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    res <- K %*% crossprod(Z, Z)
    #eig <- eigen(res, symmetric = FALSE, EISPACK = TRUE)
    eig <- magmaEigenNonsym(Xmat=res, ngpu=ngpu)
    values  <- eig$values
    vectors <- magmaQR(Xmat= (Z %*% eig$vectors), ngpu=ngpu, printInfo=FALSE)
#    vectors <- qr.Q(qr(Z %*% eig$vectors)
    return(list(values = values , vectors = vectors))
}


