 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, Zmat=NULL, numberSNPselected=0, quiet=TRUE, gamma=NULL)
 {
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, Z=Zmat, llim=-100,ulim=100)


   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

  # calculate gamma
  if(is.null(gamma)){
    print(" in here")
    # gamma has not been set
    if (length(trait) <= 700)
         gamma <- 1
    if (length(trait) > 700 & length(trait) <= 1000)
      gamma <- -0.000333 * length(trait) + 1.233333333
      # gamma <- -0.001 * length(trait) + 1.7

    if (length(trait) > 1000 & length(trait) <= 1500)
       gamma <- -0.0002 * length(trait) + 1.1

    if (length(trait) > 1500 & length(trait) <= 2000)
       gamma <- -0.0004 * length(trait) + 1.4

   if (length(trait) > 2000)
      gamma <- 0.6

  }

   extBIC <- BIC + 2 * gamma   *lchoose(geno$dim_of_ascii_M[2], numberSNPselected)  

   cat(" Gamma = ", gamma, "\n")
    return(extBIC)
 }


