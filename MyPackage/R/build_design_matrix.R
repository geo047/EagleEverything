.build_design_matrix <- function(pheno=NULL,  indxNA=NULL, fformula=NULL, quiet=TRUE  )
{
   ## internal function: use only in AM function and SummaryAM  function
   ## build design matrix given character vector fformula of column names

   ## assign model matrix X
   if(is.null(fformula))
   {  ## trait + intercept being fitted only
      if(length(indxNA) > 0){
         Xmat <- matrix(data=1, nrow=nrow(pheno[-indxNA,]), ncol=1)

      } else {
        Xmat <- matrix(data=1, nrow=nrow(pheno), ncol=1)
      }
      colnames(Xmat) <- "(Intercept)"  ## to be consistent with model.matrix
   } else {
      ## trait + fixed effects being fitted. 
     if(length(indxNA)==0 | is.null(indxNA) )
     {
        Xmat <- model.matrix(fformula, data=pheno)
     }  else {
        # there is an issue with creating Xmat when it includes
        # factors that have some of their levels removed. 
        ph <- pheno[-indxNA,]
        mat <- get_all_vars(formula=fformula, data=ph)
        for(ii in names(mat)){
           if(is.factor(ph[,ii])){
              ph[,ii] <- as.factor(as.character(ph[,ii]))
           }
        }  ## for    
        Xmat <- model.matrix(fformula, data=ph)
     } ## if else (length(indxNA)==0)
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


