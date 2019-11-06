#' @title Set the false positive rate for \code{AM}
#' @description The gamma parameter in \code{AM} controls the false positive rate of the model 
#' building process. This function uses permutation to  find the gamma value for a desired false positive rate. 
#' @param  falseposrate the desired false positive rate.
#' @param  numreps  the number of replicates upon which to base the calculation of the false 
#'                   positive rate. We have found 200 replicates to be sufficient but more is better.   
#' @param trait  the name of the column in the phenotype data file that contains the trait data. The name is case sensitive and must match exactly the column name in the phenotype data file.  This parameter must be specified. 
#' @param fformula   the right hand side formula for the fixed effects part of the model. 
#' @param geno   the R  object obtained from running \code{\link{ReadMarker}}. This must be specified. 
#' @param pheno  the R  object  obtained  from running \code{\link{ReadPheno}}. This must be specified.
#' @param map   the R object obtained from running \code{\link{ReadMap}}. If not specified, a generic map will 
#'              be assumed. 
#' @param Zmat     the R object obtained from running \code{\link{ReadZmat}}. If not specified, an identity matrix will be assumed. 
#' @param numgammas the number of equidistant gamma values from 0 to 1 for which to calculate the false positive rate of the model building process. This should not need 
#' adjusting.  
#' @param ncpu a integer  value for the number of CPU that are available for distributed computing.  The default is to determine the number of CPU automatically. 
#' @param ngpu   a integer value for the number of gpu available for computation.  The default
#'               is to assume there are no gpu available.  
#'               This option has not yet been implemented.
#' @param seed  a integer value for the starting seed for the permutations. 
#' 
#' @details

#' The false positive rate for \code{\link{AM}} is controlled by its gamma parameter. Values close to 
#' 1 (0) decreases (increases) the false positive rate of detecting SNP-trait associations. There is no 
#' analytical way of setting gamma for a specified false positive rate. So we are using permutation to do this empirically. 
#' 
#' By setting \code{falseposrate} to the desired false positive rate, this function will find the corresponding gamma value for \code{\link{AM}}. 
#' 
#' A table of other gamma values for a range of false positive rates is also given. 
#'
#' To increase the precision of the gamma estimates, increase \code{numreps}. 
#'
#'
#'
#' @seealso \code{\link{AM}}
#' @return
#' A list with the following components:
#' \describe{
#'\item{numreps:}{the number of permutations performed.}
#'\item{gamma:}{the vector of gamma values.}
#'\item{falsepos:}{the false positive rates for the gamma values.}
#'\item{setgamma:}{the gamma value that gives a false positive rate of \code{falseposrate} }
#' }
#'
#' @examples
#'   \dontrun{ 
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #-------------------------
#'   #  Example  
#'   #------------------------
#'
#'   # read the map 
#'   #~~~~~~~~~~~~~~
#'   
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'   # read marker data
#'   #~~~~~~~~~~~~~~~~~~~~
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8  gigabytes
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   # read phenotype data
#'   #~~~~~~~~~~~~~~~~~~~~~~~
#'
#'   # Read in a plain text file with data on a single trait and two covariates
#'   # The first row of the text file contains the column names y, cov1, and cov2. 
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'
#'  #  Suppose we want to perform the AM analysis at a 5% false positive rate. 
#'  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'  
#'   ans <- FPR4AM(falseposrate = 0.05,
#'                 trait = 'y',
#'                 fformula=c('cov1+cov2'),
#'                 map = map_obj,
#'                 pheno = pheno_obj,
#'                 geno = geno_obj) 
#'  
#'
#'   res <- AM(trait =  'y',
#'                 fformula=c('cov1+cov2'),
#'                 map = map_obj,
#'                 pheno = pheno_obj,
#'                 geno = geno_obj,
#'                 gamma = ans$setgamma)
#'
#'
#' }
#'
#'
FPR4AM <- function(
               falseposrate = 0.05,
               trait=trait,
               numreps = 200,
               fformula  = NULL,
               numgammas = 20,
               geno=NULL,
               pheno=NULL,
               map = NULL,
               Zmat = NULL,
               ncpu=detectCores(),
               ngpu=0,
               seed=101 
               ){

  quiet <- TRUE   ## change to FALSE if additional error checking is needed. 

  set.seed(seed)
  # need some checks in here ... 
  error.code <- check.inputs.mlam(ncpu=ncpu ,  colname.trait=trait,
                     map=map, pheno=pheno, geno=geno, Zmat=Zmat, gamma=NULL, falseposrate=falseposrate )

 if(error.code){
   message("\n The Eagle function FPR4AM has terminated with errors.\n")
   return(NULL)
 }

 ## checking if map is present. If not, generate a fake map. 
 if(is.null(map)){
   message(" Map file has not been supplied. An artificial map is being created but this map is not used in the analysis. \n")
   message(" It is only used for the reporting of results. \n")
   ## map has not been supplied. Create own map
   map <- data.frame(SNP=paste("M", 1:geno[["dim_of_ascii_M"]][2], sep=""),
                     Chr=rep(1, geno[["dim_of_ascii_M"]][2]),
                     Pos=1:geno[["dim_of_ascii_M"]][2])
  }

 indxNA_pheno <- NA
 indxNA_geno <- NA
 selected_loci <- NA
 new_selected_locus <- NA
 extBIC <- vector("numeric", 0)


 ## Turn fformula  into class formula with some checks
if(!is.null(fformula)){
 if(fformula=="")  ## added for shiny
      fformula<-NULL
 }
 if(!is.null(fformula) ){
   if(length(grep("~", fformula))==0){
      if(length(fformula)==1){
          fformula <- as.formula(paste("~", fformula, sep="") )
      }  else {
          message(" fformula has ", length(fformula), " separate terms. It should be a single statement eg. x1 + x2 + x3. \n")
          message("\n CalculateFPR  has terminated with errors.\n")
          return(NULL)
      }
   } else {
    ## problem: formula should not contain ~
    message(" Only the terms on the right hand side of the formula should be specified. \n")
    message(" Please remove the ~ from the formula. \n")
    message("\n CalculateFPR has terminated with errors.\n")
    return(NULL)
  }  ## if length grep
} ## end if(!is.null(fformula))

  ## check that terms in  formula are in pheno file
 if(!is.null(fformula)){
  res <- tryCatch(
     mat <- get_all_vars(formula=fformula, data=pheno) ,
     error = function(e)
     {
         return(TRUE)
     }
  )
  if(!is.data.frame(res))
  {
   if(res){
      message(" fformula contains terms that are not column headings in the phenotype file. \n")
      message(" Check spelling and case of terms in fformula. \n")
      message("\n  CalculateFDR has terminated with errors.\n ")
      return(NULL)
   }
  }
 }


 ## check for NA's in explanatory variables 
 ## If any, set individual's trait value to NA
 ## This means this individual will later be removed. 
 if(!is.null(fformula)){
    mat <- get_all_vars(formula=fformula, data=pheno)
    mat.of.NA  <- which(is.na(mat), arr.ind=TRUE)
  if(!is.null(dim(mat.of.NA)[1]) ){
     if(dim(mat.of.NA)[1]>0){
       pheno[unique(mat.of.NA[,1]), trait] <- NA
     }
  }
 }

 ## check for NA's in trait
 indxNA_pheno <- check.for.NA.in.trait(trait=pheno[, trait])



 ## remove missing observations from pheno  
# if(length(indxNA_pheno)>0){
#    pheno <- pheno[-indxNA_pheno, ]
#    message(" The following rows are being removed from phenotypic data  due to missing trait data: \n")
#    message(cat("             ", indxNA_pheno, "\n\n"))
# }



## If Z matrix is present, then we may not need to remove genotypes (i.e. indxNA_geno)
## check if any genotypes need to be removed. 
if(!is.null(Zmat)){
      # check for columns with 0 sums. 
      if(!any(is.na(indxNA_pheno))){
          s  <- colSums(Zmat[-indxNA_pheno,])
          indx <- which(s==0)
          if (length(indx)  > 0 )
          indxNA_geno <- indx
       }
} else {
  indxNA_geno <- indxNA_pheno
}









 # fit null model and find residuals
 # i.e Y = XB + e   calculate \hat{e}.
 pheno[, "trait"] <- pheno[, trait]  ## name change
 if(is.null(fformula)){
   ff <- as.formula("trait ~ 1")
 } else {
    # turn the formula back into a character
    fformula <- as.character(fformula)[2]  # gets rid of the ~
    ff <- paste("trait ~ ", fformula, sep=" ")
    ff <- as.formula(ff)
 }
 mod <- lm( formula=ff , data=pheno )
 res <- residuals(mod)

 if (any(is.na(indxNA_pheno))){
    pheno[, "residuals"] <- res
 } else {
    pheno[, "residuals"] <- NA
    pheno[-indxNA_pheno, "residuals"] <- res
 }

 # create big pheno: contains all permutations 
bigpheno <- matrix(data=NA, nrow=length(res) , ncol=numreps)
for(ii in 1:numreps){
  #sperm <- pheno[, "residuals"]
  sperm <- res
  sperm <- sperm[sample(1:length(sperm), length(sperm), FALSE)]
  bigpheno[, ii] <- sperm 
}
colnames(bigpheno) <- paste0("res", 1:numreps)

# bigpheno will be of the wrong dimension if there is mising 
# phenotype data (indxNA_geno)
# Need to put back these 

if (!any(is.na(indxNA_pheno))){
  # lm removes missing data by default. 
  # Therefore, bigpheno is missing rows. 
  # Adding in rows of NA's here
  id <- 1:nrow(pheno)  # size of data, including missing obs 
  id <- id[-indxNA_pheno]
  tmp <- matrix( data = NA,  nrow= nrow(pheno), ncol=ncol(bigpheno) )
  tmp[id, ] <- bigpheno
  bigpheno <- tmp
}


 ## build design matrix currentX
 currentX_null <- .build_design_matrix(pheno=bigpheno,  fformula=NULL, quiet=quiet )

 ## check currentX for solve(crossprod(X, X)) singularity
 if (any(is.na(indxNA_geno))) {
   chck <- tryCatch({ans <- solve(crossprod(currentX_null, currentX_null))},
           error = function(err){
            return(TRUE)
           })
  } else {
   chck <- tryCatch({ans <- solve(crossprod(currentX_null[-indxNA_pheno,], currentX_null[-indxNA_pheno,]))},
           error = function(err){
            return(TRUE)
           })
  }


  if(is.logical(chck)){
      if(chck){
        message(" There is a problem with the effects in fformula.\n")
        message(" These effects are causing computational instability. \n")
        message(" This can occur when there is a strong dependency between the effects.\n")
        message(" Try removing some of the effects in fformula. \n")
        message("\n  AM has terminated with errors.\n")
        return(NULL)
      }
  }



 ## calculate Ve and Vg
 Args <- list(geno=geno,
                    ncpu=ncpu,selected_loci=selected_loci,
                    quiet=quiet)
 if(!quiet)
   message(" quiet=FALSE: calculating M %*% M^t. \n")
 MMt <- do.call(.calcMMt, Args)

 if(!quiet)
   doquiet(dat=MMt, num_markers=5 , lab="M%*%M^t")
 invMMt <- chol2inv(chol(MMt))   ## doesn't use GPU
 gc()

gamma <- seq(0,1,length.out=numgammas)


# Fit null model and calculate extBIC
vc <- list()
best_ve <- rep(NA, numreps)
best_vg <- rep(NA, numreps)
MaxLike <- rep(NA, numreps)
extBIC <-   matrix(data=NA, nrow=numreps, ncol=length(gamma))

# Found that REML step gives better results, at least for small sample size

     if ( any(is.na(indxNA_pheno)) ) {

         Args <- list(y=pheno[, "residuals"] , X= currentX_null , Z=Zmat, K=MMt, ngpu=ngpu)

     }  else {
       if (is.null(Zmat)){
         # no Z 
         Args <- list(y=pheno[-indxNA_pheno , "residuals"] , X= matrix(data=currentX_null[-indxNA_pheno,], ncol=1) , Z=NULL, 
                    K=MMt[-indxNA_geno, -indxNA_geno], ngpu=ngpu)


       }  else {
         # Have to deal with Z and repeated measures
         # Case 1:  missing pheno but no missing geno
         # Case 2:  missing pheno and geno data

         # Case 1: missing pheno but no missing geno
         if (  !any(is.na(indxNA_pheno)) & any(is.na(indxNA_geno))  ){
            # check for situation where nrow(Z) is same as nrow(MMt)
            # this causes issues for EMMA later
            if (  !any(colSums(Zmat[-indxNA_pheno, ] )> 1) ){
                Args  <- list(y=pheno[-indxNA_pheno , "residuals"] , X= matrix(data=currentX_null[-indxNA_pheno,], ncol=1) , Z=NULL, 
                              K=MMt[-indxNA_geno, -indxNA_geno], ngpu=ngpu)


            } else {
              Args <- list(y=pheno[-indxNA_pheno , "residuals"] , X= matrix(data=currentX_null[-indxNA_pheno,], ncol=1) , Z=Zmat[-indxNA_pheno, ],
                           K=MMt[-indxNA_geno, -indxNA_geno], ngpu=ngpu)

            }  ## end if nrow(Zmat) == nrow(MMt)

          }   ## end if Case 1

          # Case 2:  missing pheno and missing geno 
          if (  !any(colSums(Zmat[-indxNA_pheno, -indxNA_geno] )> 1) ){
             # this captures case of Z turning into an identify matrix
             Args  <- list(y=pheno[-indxNA_pheno , "residuals"] , X= matrix(data=currentX_null[-indxNA_pheno,], ncol=1) , Z=NULL ,
                           K=MMt[-indxNA_geno, -indxNA_geno], ngpu=ngpu)



          } else {
             Args  <- list(y=pheno[-indxNA_pheno , "residuals"] , X= matrix(data=currentX_null[-indxNA_pheno,], ncol=1)  , 
                           Z=Zmat[-indxNA_pheno, -indxNA_geno]  , K=MMt[-indxNA_geno, -indxNA_geno], 
                           ngpu=ngpu)
         } ## end if  Case 2




       }  ## end  if (is.null(Zmat))


     } # end  if ( any(is.na(indxNA_pheno)) )


   res_full  <- do.call(emma.REMLE, Args)
   vc <- list("vg"=res_full$vg, "ve"=res_full$ve   )

rep(NA, numreps)

 for (ii in 1:numreps){
   #vc[[ii]] <- .calcVC(trait=bigpheno[, ii], Zmat=Zmat, currentX=currentX_null,MMt=MMt, ngpu=ngpu)
   #gc()
   #best_ve[ii] <- vc[[ii]][["ve"]]
   #best_vg[ii] <- vc[[ii]][["vg"]]
   #MaxLike[ii] <- vc[[ii]][["ML"]]



   # Approximation step here - vc being computed once on the original unpermuted trait (residuals). 
   # Exact code is what is commented out above where vc computed on each permutation of residuals.  
   best_ve[ii] <- vc[["ve"]]
   best_vg[ii] <- vc[["vg"]]
   #MaxLike[ii] <- vc[["ML"]]
   gc()
 } ## for ii


# Since there are no fixed effects and we have permuted Y, then 
# we only need to do this once, for a single rep and all the rest
# will have the same null extBIC value. 
 if (  any(is.na(indxNA_pheno)) ){
 Args <- list("trait"= bigpheno[,1], "currentX"=currentX_null, "geno"=geno, "MMt"=MMt,
                       "Zmat"=Zmat, "numberSNPselected"=0, "quiet"=quiet, "gamma"=gamma)
 }
  if ( !any(is.na(indxNA_pheno))  &  any(is.na(indxNA_geno)) ){
    if ( all(colSums( Zmat[-indxNA_pheno,]) == 1 ) & all(rowSums( Zmat[-indxNA_pheno,]) == 1 ) ){
          Zmatrix = NULL
    } else {
      Zmatrix = Zmat[-indxNA_pheno,]
    }
    Args <- list("trait"= bigpheno[-indxNA_pheno ,1], "currentX"=matrix(data=currentX_null[-indxNA_pheno,], ncol=1)  , "geno"=geno, "MMt"=MMt,
                       "Zmat"=Zmatrix , "numberSNPselected"=0, "quiet"=quiet, "gamma"=gamma)
  }

 if ( !any(is.na(indxNA_pheno))  & !any(is.na(indxNA_geno))  ){
    if ( all(colSums( Zmat[-indxNA_pheno, -indxNA_geno]) == 1 ) & all(rowSums( Zmat[-indxNA_pheno, -indxNA_geno]) == 1 ) ){
          Zmatrix = NULL
    } else {
      Zmatrix = Zmat[-indxNA_pheno, -indxNA_geno]
    }
    Args <- list("trait"= bigpheno[-indxNA_pheno ,1], "currentX"=matrix(data=currentX_null[-indxNA_pheno,], ncol=1)  , "geno"=geno, "MMt"=MMt[-indxNA_geno, -indxNA_geno],
                       "Zmat"=Zmatrix, "numberSNPselected"=0, "quiet"=quiet, "gamma"=gamma)
 }


 extBIC <-     do.call(.calc_extBIC_MLE, Args)
     gc()

extBIC <- matrix(data=extBIC, nrow=numreps, ncol=length(gamma)) # formed matrix of null extBIC values


 extBIC_alternate <- matrix(data=NA, nrow=numreps, ncol=length(gamma))
 ## => only need to do once. 
 ## Looks at the stability of the MMt calculation especially if there are near identical rows of data in M
 error_checking <- FALSE
 if (!quiet ) error_checking <- TRUE

 MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking,
                                                              ngpu=ngpu , message=message)

      if (  any(is.na(indxNA_pheno)) ){
         H <- calculateH(MMt=MMt, varE=best_ve[ii], varG=best_vg[ii], Zmat=Zmat, message=message )
         P <- calculateP(H=H, X=currentX_null , message=message)
      hat_a <- calculate_reduced_a_batch(Zmat=Zmat, varG=best_vg[ii], P=P,
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]],
                       y=bigpheno , quiet = quiet , message=message)
      var_hat_a    <- calculate_reduced_vara(Zmat=Zmat, X=currentX_null, 
                                             varE=best_ve[ii], varG=best_vg[ii], invMMt=invMMt,
                                             MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]],
                                             quiet = quiet, message=message )
      } else {

         H <- calculateH(MMt=MMt, varE=best_ve[ii], varG=best_vg[ii], Zmat=Zmat[-indxNA_pheno, ], message=message )
         P <- calculateP(H=H, X=matrix(data=currentX_null[-indxNA_pheno, ], ncol=1) , message=message)
      hat_a <- calculate_reduced_a_batch(Zmat=Zmat[-indxNA_pheno, ] , varG=best_vg[ii], P=P,
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]],
                       y=bigpheno[-indxNA_pheno,  ]  , quiet = quiet , message=message)

      var_hat_a    <- calculate_reduced_vara(Zmat=Zmat[-indxNA_pheno, ], X=matrix(data=currentX_null[-indxNA_pheno, ], ncol=1) , 
                                             varE=best_ve[ii], varG=best_vg[ii], invMMt=invMMt,
                                             MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]],
                                             quiet = quiet, message=message )

      }
 

   
      a_and_vara  <- calculate_a_and_vara_batch(numreps = numreps, 
                                          geno = geno,
                                          selectedloci = selected_loci,
                                          invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]],
                                          transformed_a=hat_a ,
                                          transformed_vara=var_hat_a,
                                          quiet=quiet, message=message,  indxNA_geno = indxNA_geno)




for(ii in 1:numreps){
       if(ii %% 4 == 0 )
          message("-", appendLF=FALSE)
       if(ii %% 4 == 1 )
          message("\\", appendLF=FALSE)
       if(ii %% 4 == 2 )
          message("|", appendLF=FALSE)
       if(ii %% 4 == 3 )
          message("/", appendLF=FALSE)
 


      # Select new locus : find_qtl function but with calculateMMt_sqrt_and_sqrtinv
      # moved outside the function for computational gain with permuted samples
     

      ## outlier test statistic
###   tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
      tsq <- a_and_vara[["a"]][, ii]**2/a_and_vara[["vara"]]

      indx <- which(tsq == max(tsq, na.rm=TRUE))   ## index of largest test statistic. 
                                                   ## However, need to account for other loci 
                                                   ## already having been removed from M which 
                                                   ## affects the indexing

       ## taking first found qtl
       midpoint <- 1
       if (length(indx)>2){
          midpoint <- trunc(length(indx)/2)+1
       }
       indx <- indx[midpoint]

       new_selected_locus  <- seq(1, geno[["dim_of_ascii_M"]][2])  ## 1:ncols
       new_selected_locus <- new_selected_locus[indx]

   
       # Fit alternate model
       currentX <- currentX_null # initialising to base model  
#       if (  any(is.na(indxNA_pheno)) ){
       currentX <- constructX(Zmat=Zmat, fnameMt=geno[["asciifileMt"]], currentX=currentX, 
                              loci_indx=new_selected_locus,
                              dim_of_Mt=geno[["dim_of_ascii_Mt"]],
                              map=map )
#       } else {
#       currentX <- constructX(Zmat=Zmat , fnameMt=geno[["asciifileMt"]], currentX=currentX , 
#                              loci_indx=new_selected_locus,
#                              dim_of_Mt=geno[["dim_of_ascii_Mt"]],
#                              map=map )
#
#      }


     if (  any(is.na(indxNA_pheno)) ){
            Args <- list("trait"= bigpheno[,ii], "currentX"=currentX, "geno"=geno, "MMt"=MMt,
                       "Zmat"=Zmat, "numberSNPselected"=1, "quiet"=quiet, "gamma"=gamma)
     }
     if ( !any(is.na(indxNA_pheno))  &  any(is.na(indxNA_geno)) ){
    if ( all(colSums( Zmat[-indxNA_pheno,]) == 1 ) & all(rowSums( Zmat[-indxNA_pheno,]) == 1 ) ){
          Zmatrix = NULL
    } else {
      Zmatrix = Zmat[-indxNA_pheno,]
    }

             Args <- list("trait"= bigpheno[-indxNA_pheno ,ii], "currentX"=currentX[-indxNA_pheno,]  , "geno"=geno, "MMt"=MMt,
                       "Zmat"=Zmatrix, "numberSNPselected"=1, "quiet"=quiet, "gamma"=gamma)
      }

     if ( !any(is.na(indxNA_pheno))  & !any(is.na(indxNA_geno))  ){
    if ( all(colSums( Zmat[-indxNA_pheno,-indxNA_geno]) == 1 ) & all(rowSums( Zmat[-indxNA_pheno,-indxNA_geno]) == 1 ) ){
          Zmatrix = NULL
    } else {
      Zmatrix = Zmat[-indxNA_pheno,-indxNA_geno]
    }

             Args <- list("trait"= bigpheno[-indxNA_pheno ,ii], "currentX"=currentX[-indxNA_pheno,]  , "geno"=geno, "MMt"=MMt[-indxNA_geno, -indxNA_geno],
                       "Zmat"=Zmatrix, "numberSNPselected"=1, "quiet"=quiet, "gamma"=gamma)
      }





 extBIC_alternate[ii, ] <-     do.call(.calc_extBIC_MLE, Args)

}



# calculate FDR for gamma value
mat <- (extBIC_alternate < extBIC)
falsepos <- colSums(mat)/numreps

cat("\n\n")
cat(" Table: Empirical false positive rates, given gamma value for model selection. \n\n")
cat("  Gamma    |  False Pos Rate  \n")
cat(" ---------------------------- \n")
for(ii in 1:length(gamma)){
cat( gamma[ii], " | ", round(falsepos[ii], 3), "\n")
} 
cat(" ----------------------------- \n")

# find best gamma value 

d <- abs(falsepos - falseposrate)
indx <- which(min(d) == d)
if(length(indx) > 1){
    indx <- max(indx)   ## picking most conservative value if there are multiple values
}
 setgamma <- gamma[indx]

cat(" For a false positive rate of ", falseposrate, " set the gamma parameter in the AM function to ", setgamma, "\n")


 return(list(numreps=numreps, gamma=gamma, falsepos=falsepos, setgamma = setgamma    ))



} ## end function FPR4AMnew





