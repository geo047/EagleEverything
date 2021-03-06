  #.find_qtl <- function(Zmat=NULL, geno, availmemGb,  selected_loci, MMt, invMMt, best_ve, best_vg,
  #                     currentX,  ncpu, quiet, trait, ngpu  )
  .find_qtl <- function(Zmat=NULL, geno, availmemGb,  selected_loci, MMt, invMMt, best_ve, best_vg,
                       currentX,  ncpu, quiet, trait, ngpu, itnum )
  {
    ##  internal function: use by   AM
    print(" inside find.qtl")
    H <- calculateH(MMt=MMt, varE=best_ve, varG=best_vg, Zmat=Zmat, message=message )
    print("H")
    print(H)
    if(!quiet)
        doquiet(dat=H, num_markers=5, lab="H")
    print(" currentX")
    print(currentX)

    P <- calculateP(H=H, X=currentX , message=message)
    if(!quiet)
        doquiet(dat=P, num_markers=5 , lab="P")
    print("P")
    print(P)
    rm(H)
    gc()




    ## Looks at the stability of the MMt calculation especially if there are near identical rows of data in M
    error_checking <- FALSE
    if (!quiet )
       error_checking <- TRUE
    MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking,
                              ngpu=ngpu , message=message)
    print("MMt_sqrt_and_sqrtinv[[sqrt_MMt]]")
    print(MMt_sqrt_and_sqrtinv[["sqrt_MMt"]])
   

    if(!quiet){
       doquiet(dat=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], num_markers=5, lab="sqrt(M %*% M^t)")
       doquiet(dat=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]], num_markers=5, lab="sqrt(M %*% M^t)^-1")
    }
    if(!quiet ){
      message(" quiet =", quiet, ": beginning calculation of the BLUP estimates for dimension reduced model. \n")
    }
    hat_a <- calculate_reduced_a(Zmat=Zmat, varG=best_vg, P=P,
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]],
                       y=trait, quiet = quiet , message=message)
    if(!quiet)
       doquiet(dat=hat_a, num_markers=5, lab="BLUPs")
    print("hat_a")
    print(hat_a)


     rm(P)
     gc()

    if(!quiet ){
      message(" quiet = ", quiet, ": beginning calculation of the standard errors  of BLUP estimates for dimension reduced model. \n")
    }

    var_hat_a    <- calculate_reduced_vara(Zmat=Zmat, X=currentX, varE=best_ve, varG=best_vg, invMMt=invMMt,
                                                MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]],
                                                quiet = quiet, message=message )
    if(!quiet)
             doquiet(dat=var_hat_a, num_markers=5, lab="SE of BLUPs")



     gc()
    if(!quiet ){
      message(" quiet = ", quiet, ": beginning calculation of BLUPS and their standard errors for full model. \n")
    }

     a_and_vara  <- calculate_a_and_vara(geno = geno,
                                         maxmemGb=availmemGb,
                                            selectedloci = selected_loci,
                                            invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]],
                                            transformed_a=hat_a,
                                            transformed_vara=var_hat_a,
                                            quiet=quiet, message=message)
     if(!quiet){
        doquiet(dat=a_and_vara[["a"]], num_markers=5, lab="BLUPs for full model")
        doquiet(dat=a_and_vara[["vara"]], num_markers=5, lab="SE of BLUPs for full model")
     }

    ## outlier test statistic
    if (!quiet )
        message(" quiet = ", quiet, ": beginning calculation of outlier test statistics. \n")
    print(" checking tsq ")
    print(a_and_vara[["a"]][1:5])
    print(a_and_vara[["vara"]][1:5])
    tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
    if(!quiet)
       doquiet(dat=tsq, num_markers=5, lab="outlier test statistic")


   # saving test statistic to disc so that we can plot these results
#   if(.Platform$OS.type == "unix") {
#       tmpfile <- paste(tempdir(), "/", paste("tsq", itnum , ".RData", sep="") , sep="")
#       save(tsq, file=tmpfile)
#     } else {
#       tmpfile <- paste(tempdir()  , "\\", paste("tsq", itnum , ".RData", sep="")     , sep="")
#       save(tsq, file=tmpfile)
#     }

 ## AWG tmp
# tmpfile <- paste("/home/geo047/"   , "/", paste("tsq", itnum , ".RData", sep="")     , sep="")
#  save(tsq, file=tmpfile)





    indx <- which(tsq == max(tsq, na.rm=TRUE))   ## index of largest test statistic. However, need to account for other loci 
                                         ## already having been removed from M which affects the indexing

    ## taking first found qtl
    midpoint <- 1
    if (length(indx)>2){
      midpoint <- trunc(length(indx)/2)+1
    } 
    indx <- indx[midpoint]

    orig_indx <- seq(1, geno[["dim_of_ascii_M"]][2])  ## 1:ncols
    res <- list()
    res[["orig_indx"]] <- orig_indx[indx]
    res[["outlierstat"]] <- tsq
    return(res)
    #return(orig_indx[indx])
}



