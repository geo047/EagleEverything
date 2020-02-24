 .calc_extBIC <- function( ML=0 , trait=NULL, currentX=NULL,   geno=NULL,  numberSNPselected=0, quiet=TRUE, gamma=NULL)
 {
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only

   #res_p <- emma.MLE(y=trait, X= currentX , K=MMt, Z=Zmat, llim=-100,ulim=100)
   #BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

   BIC <- -2 * ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

  # calculate gamma
  if(is.null(gamma)){
     gamma <- 1
     # found this to be anti-conservative when sample size is small
     #lambda <- log(geno$dim_of_M[2])/log(length(trait))
     #gamma <- 1-(1/(2*lambda))
   } ## outer if

  if(length(gamma)==1){
   extBIC <- BIC + 2 * gamma   *lchoose(geno$dim_of_M[2], numberSNPselected)  
   if(!quiet){
      cat(" Gamma = ", gamma, "\n")
   }
   } else {
     extBIC <- rep(NA, length(gamma))
     for(ii in 1:length(gamma)){
        extBIC[ii] <- BIC + 2 * gamma[ii]  *lchoose(geno$dim_of_M[2], numberSNPselected)
     }

   }


    return(extBIC)
 }


