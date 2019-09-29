emma.eigen.L.wo.Z <- function (K, ngpu=1)
{
   # eig <- eigen(K, symmetric = TRUE)
    eig <- magmaEigen(Xmat=K, ngpu=ngpu, printInfo=TRUE)
    return(list(values = eig$values, vectors = eig$vectors))
}


