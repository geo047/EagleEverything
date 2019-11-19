 emma.eigen.L.w.Z <- function (Z, K, complete = TRUE )
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    res <- K %*% crossprod(Z, Z)

    eig <- eigen(res, symmetric = FALSE, EISPACK = TRUE)

    values  <- eig$values

   if( nrow( (Z %*% eig$vectors) ) != ncol ( (Z %*% eig$vectors) )){
    vectors <- qr.Q(qr(Z %*% eig$vectors) )
   } else {
    vectors <- magmaQR(Xmat= (Z %*% eig$vectors),  printInfo=FALSE)
   }
    return(list(values = values , vectors = vectors))
}


