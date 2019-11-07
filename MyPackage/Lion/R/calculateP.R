
calculateP  <- function(H=NULL, X=NULL, ngpu=1)
{
  ## internal function to AM
  ## R function to calculate P matrix
  ## Args:
  ##       H is the variance matrix
  ##       X is the design matrix supplied by the user
  ##       ngpu for the number of gpu
  ## Returns:
  ##   matrix object P

  if(is.null(H)){
    message(" H must be specified.")
    return(NULL)
  }
  if(is.null(X)){
    message(" A design matrix has not be specified. ")
    return(NULL)
  }

   if(nrow(H) != nrow(X)){
      message(" The number of rows in H and X are not the same.")
    return(NULL)
  }

  Hinv <- chol2inv(chol(H))
 #P <- Hinv - Hinv %*% X %*% magmaSolve( Xmat= (t(X) %*% Hinv %*% X ), ngpu=ngpu, printInfo=FALSE)  %*% t(X) %*% Hinv
 P <- Hinv - Hinv %*% X %*% chol2inv(chol( t(X) %*% Hinv %*% X ))     %*% t(X) %*% Hinv

  return(P)

}


