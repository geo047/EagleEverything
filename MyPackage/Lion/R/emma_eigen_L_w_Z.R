 emma.eigen.L.w.Z <- function (Z, K, complete = TRUE )
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    res <- K %*% crossprod(Z, Z)

   doMagmaEigen <- FALSE
   if (computer$ngpu == 1) {
     if ( nrow(res) > 3000 )
           doMagmaEigen <- TRUE
   } else if (computer$ngpu == 2) {
     if (nrow(res) > 3400 )
           doMagmaEigen <- TRUE
   } else if (computer$ngpu == 3) {
     if (nrow(res) > 4100 )
           doMagmaEigen <- TRUE
  } else if (computer$ngpu > 3) {
     if (nrow(res) > 4100)
           doMagmaEigen <- TRUE
  } else {
    message(" Error: the number of gpu needs to be set to a sensible number. \n")
    return(NULL)
 }
   if (doMagmaEigen){
    eig <- magmaEigenNonsym(Xmat=res)
   } else {
    eig <- eigen(res, symmetric = FALSE, EISPACK = TRUE)
   }

    values  <- eig$values

   if( nrow( (Z %*% eig$vectors) ) != ncol ( (Z %*% eig$vectors) )){
    vectors <- qr.Q(qr(Z %*% eig$vectors) )
   } else {
    vectors <- magmaQR(Xmat= (Z %*% eig$vectors),  printInfo=FALSE)
   }
    return(list(values = values , vectors = vectors))
}


