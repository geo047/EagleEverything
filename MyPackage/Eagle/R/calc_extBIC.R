 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, Zmat=NULL, numberSNPselected=0, quiet=TRUE, gamma=NULL)
 {
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, Z=Zmat, llim=-100,ulim=100)


   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

  # calculate gamma
  if(is.null(gamma)){
     lambda <- log(geno$dim_of_ascii_M[2])/log(length(trait))
     gamma <- 1-(1/(2*lambda))
   }


   extBIC <- BIC + 2 * gamma   *lchoose(geno$dim_of_ascii_M[2], numberSNPselected)  
   cat(" Gamma = ", gamma, "\n")
    return(extBIC)
 }


