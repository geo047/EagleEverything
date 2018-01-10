 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, Zmat=NULL, numberSNPselected=0, quiet=TRUE)
 {
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, Z=Zmat, llim=-100,ulim=100)


   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

   ## calculate lambda for extended BIC
   #lambda <- log(geno$dim_of_ascii_M[2])/log(length(trait))
   #gamma <- 1-(1/(2*lambda))
   #extBIC <- BIC + 2 * gamma *lchoose(geno$dim_of_ascii_M[2], numberSNPselected)  # anti-conservative
   #extBIC <- BIC + 2 *  1.0   *lchoose(geno$dim_of_ascii_M[2], numberSNPselected)  
   extBIC <- BIC + 2 *  (1 - ( log(geno$dim_of_ascii_M[1])/(4 * log(geno$dim_of_ascii_M[2])) ) )   *lchoose(geno$dim_of_ascii_M[2], numberSNPselected)  
    return(extBIC)
 }


