  .calcVC <- function(trait, Zmat, currentX, MMt, eigenMMt, ngpu)
  {
    ## internal function: used by AM and SummaryAM
    ## perform likelihood ratio test for variance component Var_g
    # res_full <- emma.REMLE(y=trait, X= currentX , Z=Zmat, K=MMt, ngpu=ngpu)


 
    res_full <- emma.MLE(y=trait, X= currentX , Z=Zmat, K=MMt, ngpu=ngpu )
    return(list("vg"=res_full$vg, "ve"=res_full$ve, "ML"=res_full$ML   ))

  }


