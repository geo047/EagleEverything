 emma.eigen.L.w.Z <- function (Z, K, complete = TRUE, ngpu=1)
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    res <- K %*% crossprod(Z, Z)

   doMagmaEigen <- FALSE
   if (ngpu == 1) {
     if ( nrow(res) > 3000 )
           doMagmaEigen <- TRUE
   } else if (ngpu == 2) {
     if (nrow(res) > 3400 )
           doMagmaEigen <- TRUE
   } else if (ngpu == 3) {
     if (nrow(res) > 4100 )
           doMagmaEigen <- TRUE
  } else if (ngpu > 3) {
     if (nrow(res) > 4100)
           doMagmaEigen <- TRUE
  } else {
    message(" Error: the number of gpu needs to be set to a sensible number. \n")
    return(NULL)
 }

   if (doMagmaEigen){
    eig <- magmaEigenNonsym(Xmat=res, ngpu=ngpu)
   } else {
    eig <- eigen(res, symmetric = FALSE, EISPACK = TRUE)
   }

    values  <- eig$values
    vectors <- magmaQR(Xmat= (Z %*% eig$vectors), ngpu=ngpu, printInfo=FALSE)
#    vectors <- qr.Q(qr(Z %*% eig$vectors)
    return(list(values = values , vectors = vectors))
}


