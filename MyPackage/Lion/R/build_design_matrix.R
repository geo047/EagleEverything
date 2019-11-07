.build_design_matrix <- function(pheno=NULL,   fformula=NULL, quiet=TRUE  )
{
   ## internal function: use only in AM function and SummaryAM  function
   ## build design matrix given character vector fformula of column names

   ## assign model matrix X
   if(is.null(fformula))
   {  ## trait + intercept being fitted only
      Xmat <- matrix(data=1, nrow=nrow(pheno), ncol=1)
      colnames(Xmat) <- "(Intercept)"  ## to be consistent with model.matrix
   } else {
      ## trait + fixed effects being fitted. 
      current.na.action <- options('na.action')[[1]]
      options(na.action='na.pass')   # this leaves NA's in the model matrix, preserving its correct dim
      Xmat <- model.matrix(fformula, data=pheno)
      options(na.action=current.na.action)
   }

 if (!quiet ){
   message("Dimension of design matrix, before addition of marker fixed effects is ", nrow(Xmat), " rows and ", ncol(Xmat), " columns.\n")
 }
if(!is.matrix(Xmat))
   Xmat <- matrix(data=Xmat, ncol=1)

## check that matrix doesn't all contain the same value
indx <- NULL
if (ncol(Xmat) > 1 ){
  for(ii in 2:ncol(Xmat)){
    # first column has intercept
    u <- length(unique(Xmat[, ii]))
    if(u == 1){
      indx <- c(indx, ii)
    }

  }
  if(length(indx) > 0)
     Xmat <- Xmat[, -indx]
}  







  return(Xmat)
}


