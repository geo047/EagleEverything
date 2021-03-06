#' @title multiple-locus Association Mapping 
#' @description \code{AM} performs  association mapping within a multiple-locus linear mixed model framework. 
#' \code{AM}  finds the best set of 
#' marker loci in strongest association with a trait while simultaneously accounting for any fixed effects and the genetic background.     
#' @param trait  the name of the column in the phenotype data file that contains the trait data. The name is case sensitive and must match exactly the column name in the phenotype data file. 
#' @param fformula   the right hand side formula for the fixed effects.   See below for details. 
#'                        If
#'                        not specified, only an overall mean will be fitted.
#' @param geno   the R  object obtained from running \code{\link{ReadMarker}}. This must be specified. 
#' @param pheno  the R  object  obtained  from running \code{\link{ReadPheno}}. This must be specified.
#' @param map   the R object obtained from running \code{\link{ReadMap}}. If not specified, a generic map will 
#'              be assumed. 
#' @param Zmat     the R object obtained from running \code{\link{ReadZmat}}. If not specified, an identity matrix will be assumed. 
#' @param ncpu a integer  value for the number of CPU that are available for distributed computing.  The default is to determine the number of CPU automatically. 
#' @param ngpu   a integer value for the number of gpu available for computation.  The default
#'               is to assume there are no gpu available.  This option has not yet been implemented.
#' @param  quiet      a logical value. If set to \code{FALSE}, additional runtime output is printed. 
#' This is useful for error checking and monitoring the progress of a large analysis. 
#' @param maxit     an integer value for the maximum number of forward steps to be performed.  This will rarely need adjusting. 
#' @param fixit     a boolean value. If TRUE, then \code{maxit} iterations are performed, regardless of the value of the model fit value extBIC. If FALSE, 
#' then the model building process is stopped when extBIC increases in value. 
#' @param lambda     a value between 0 and 1 for the regularization parameter for the extBIC. Values close to 0 lead to an anti-conservative test. Values close to 1 lead to a  
#' more conservative test. If this value is left unspecified, a default value of 1 is assumed. See \code{\link{FPR4AM}} for an empirical approach for setting the  lambda value. 
#'
#' @details
#'
#' This function is used to perform genome-wide association mapping. The phenotypic and SNP data should already be read in prior to running this function 
#' (see below for examples).  \code{AM} builds the linear mixed model iteratively, via forward selection. It is through this model building process that we 
#' identify the SNP-trait associations.  We use the extended BIC (extBIC) to decide on the 'best' model and when to stop looking for a better model. 
#' The conservativeness of extBIC can be adjusted.  If the \code{lambda} parameter is left at is default setting, then \code{AM} is run in its most 
#' conservative state (i.e. false positives are minimized but this also decreases the chance of true positives). 
#'
#' When interested in running  \code{AM}  at a certain false positive rate, use \code{\link{FPR4AM}}. This function uses permutation to 
#' find the lambda value for a desired false positive rate for \code{AM}. 
#'
#' Below are some examples of how to use \code{AM} for genome-wide association mapping of data. 
#'
#' \subsection{How to perform a basic AM analysis}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data are contained in the file geno.txt which is a plain space separated
#' text file with no column headings. The file is located in the current working directory. 
#' It contains numeric genotype values 0, 1, and 2 for snp genotypes
#' AA, AB, and BB, respectively. It also contains the  value X for a missing genotype. }
#' \item{}{the phenotype data is contained in the file pheno.txt which is a plain space
#' separated text file containing a single column with the trait data. The first row of the file 
#' has the column heading 'y'. This file does not contain any missing data.  
#' The file is located in the current working directory.}
#' \item{}{there is no map data.}
#' }
#'
#'  To analyse these data, we would use the following three functions (the parameters can be specified in any order, as well as the 
#' functions as long as AM is run last):
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename='geno.txt', AA=0, AB=1, BB=2, type="text", missing='X')
#'   
#'   pheno_obj <- ReadPheno(filename='pheno.txt')
#'
#'  # since lambda is not specified, this will run AM conservatively (where the false positive rate is lowest).
#'   res <- AM(trait='y', geno=geno_obj, pheno=pheno_obj)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{How to perform a more complicated AM analysis where the false positive rate is 5\%}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data are contained in the file geno.ped which is a 'PLINK' ped file. See
#' \code{\link{ReadMarker}} for details. The file is located in /my/dir. Let's assume 
#' the file is large, say 50 gigabytes,   and our computer only has 32 gigabytes of RAM.}
#' \item{}{the phenotype data is contained in the file pheno.txt which is a plain space
#' separated text file with  six columns. The first row of the file contains the column headings. 
#' The first column is a trait and is labeled y1.
#' The second column is another trait and is labeled y2. The third and fourth columns 
#' are nuisance variables and are labeled cov1 and cov2. The fifth and sixth columns
#' are the first two principal components to account for population substructure and are 
#' labeled pc1 and pc2. The file contains missing data that are coded as 99. 
#' The file is located in /my/dir.}
#' \item{}{the map data is contained in the file map.txt, is also located in 
#'  /my/dir, and the first row has the column headings.}
#' \item{}{An 'AM' analysis is performed where the trait of interest is y2, and
#' the fixed effects part of the model is cov1 + cov2 + pc1 + pc2, }
#' } 
#'
#'  To analyse these data, we would run the following:
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename='/my/dir/geno.ped', type='PLINK', availmemGb=32)
#'   
#'   pheno_obj <- ReadPheno(filename='/my/dir/pheno.txt', missing=99)
#'
#'   map_obj   <- ReadMap(filename='/my/dir/map.txt')
#'
#'   # FPR4AM calculates the lambda value corresponding to a desired false positive rate of 5\%
#'   ans <- FPR4AM(falseposrate=0.05, numreps=200, trait='y2', fformula=c('cov1 + cov2 + pc1 + pc2'), 
#'             geno=geno_obj, pheno=pheno_obj, map=map_obj)
#'
#'   # performs association mapping with a 5\% false positive rate
#'   res <- AM(trait='y2', fformula=c('cov1 + cov2 + pc1 + pc2'), 
#'             geno=geno_obj, pheno=pheno_obj, map=map_obj,  lambda=ans$setlambda)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{How to perform an analysis where individuals have multiple observations}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data are contained in the file geno.ped which is a 'PLINK' ped file. See
#' \code{\link{ReadMarker}} for details. The file is located in /my/dir. Let's assume 
#' the file is large, say 50 gigabytes,   and our computer only has 32 gigabytes of RAM.}
#' \item{}{the phenotype data is contained in the file pheno.txt which is a plain space
#' separated text file with  six columns. The first row of the file contains the column headings. 
#' The first column is a trait and is labeled y1.
#' The second column is another trait and is labeled y2. The third and fourth columns 
#' are nuisance variables and are labeled cov1 and cov2. The fifth and sixth columns
#' are the first two principal components to account for population substructure and are 
#' labeled pc1 and pc2. The file contains missing data that are coded as 99. 
#' The file is located in /my/dir.}
#' \item{}{the Z matrix data are contained in the file Zmatrix.txt. The file is located in /my/dir.This file is
#' a design matrix that only contains zeros and ones where each row must contain only a single one in the column that matches 
#' the individual's trait value to their corresponding genotype.}
#' \item{}{the map data is contained in the file map.txt, is also located in 
#'  /my/dir, and the first row has the column headings.}
#' \item{}{An 'AM' analysis is performed where the trait of interest is y2,  and
#' the fixed effects part of the model is cov1 + cov2 + pc1 + pc2.}
#' } 
#'
#'  To analyse these data, we would run the following:
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename='/my/dir/geno.ped', type='PLINK', availmemGb=32)
#'   
#'   pheno_obj <- ReadPheno(filename='/my/dir/pheno.txt', missing=99)
#'
#'   map_obj   <- ReadMap(filename='/my/dir/map.txt')
#'
#'   Zmat_obj  <- ReadZmat(filename='/my/dir/Zmatrix.txt')
#'
#'   res <- AM(trait='y2', fformula=c('cov1 + cov2 + pc1 + pc2'), 
#'             geno=geno_obj, pheno=pheno_obj, map=map_obj, Zmat=Zmat_obj )
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{Dealing with missing marker data}{
#'
#' \code{AM} can tolerate some missing marker data. However, ideally, 
#' a specialized genotype imputation program such as  'BEAGLE'  or 'PHASE2', should be 
#' used to impute the missing marker data before being read into 'Eagle'.  
#'
#' }
#'
#' \subsection{Dealing with missing trait data}{
#'
#'  \code{AM} deals automatically with individuals with missing trait data. 
#' These individuals are removed  from the analysis and a warning message is generated.
#' }
#' 
#' \subsection{Dealing with missing explanatory variable values}{
#'
#' \code{AM} deals automatically with individuals with missing explanatory variable values. 
#' These individuals are removed from the analysis and a warning message is generated
#' }
#'
#' \subsection{Error Checking}{
#'
#' Most errors occur when reading in the data. However, as an extra precaution, if \code{quiet=FALSE}, then additional 
#' output is printed during the running of \code{AM}. If \code{AM} is failing, then this output can be useful for diagnosing 
#' the problem. 
#'}
#'
#'
#'
#'
#' @seealso \code{\link{FPR4AM}} , \code{\link{ReadMarker}}, \code{\link{ReadPheno}},  \code{\link{ReadZmat}}, and \code{\link{ReadMap}}
#'
#' @return
#' A list with the following components:
#' \describe{
#'\item{trait:}{column name of the trait being used by 'AM'.}
#'\item{fformula:}{the fixed effects part of the linear mixed model.}
#'\item{indxNA_pheno:}{a vector containing the row indexes of the phenotyic data  that have been removed from the analysis.}
#'\item{indxNA_geno:}{a vector containing the row indexes of those genotypes that have been removed from the analysis 
#' due to missing data. }
#' \item{Mrk:}{a vector with the names of the snp in strongest and significant association with the trait.If no loci are found to be 
#' significant, then this component is \code{NA}.}
#' \item{Chr:}{the chromosomes on which the identified snp lie.}
#' \item{Pos:}{the map positions for the identified snp.}
#' \item{Indx:}{the column indexes in the marker file of the identified snp.} 
#' \item{ncpu:}{number of cpu used for the calculations.}
#' \item{availmemGb:}{amount of RAM in gigabytes that has been set by the user.}
#' \item{quiet:}{ boolean value of the parameter.}
#' \item{extBIC:}{numeric vector with the extended BIC values for the loci  found to be in  significant association with the trait.}
#' \item{lambda}{the numeric value of the parameter.}
#'}
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
#'  # Performing multiple-locus genome-wide association mapping with a model 
#'  #    with fixed effects cov1 and cov2 and an intercept. The intercept 
#'  #    need not be specified as it is assumed. 
#'  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'  
#'   res <- AM(trait = 'y',
#'                            fformula=c('cov1+cov2'),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj )
#' }
#'
AM <- function(trait=NULL, 
               fformula  = NULL,
               geno=NULL, 
               pheno=NULL, 
               map = NULL,
               Zmat = NULL,
               ncpu=detectCores(),
               ngpu=0,
               quiet=TRUE,
               maxit=40,
               fixit=FALSE,
               lambda=1 
               ){

 ## Core function for performing whole genome association mapping with EMMA
 ## Args
 ## ncpu        number of cores available for computation
 ## pheno           data frame 
 ##                 remaining columns are explanatory variables to include in the model. If a numeric vector, then it 
 ##                 is only a response to be fitted. 
 ## geno            if geno is a matrix or data frame, then the user has not ReadMarker and a bin packed file
 ##                 has not been created. If it is a character string, then it is the file location of the binary packed files. 
 ## maxit           maximum number of qtl to include in the model
 ## ngpu            number of gpu available for computation

 assign("ngpu", 0 , envir=computer)


 if(is.null(lambda))
    lambda <- 1


 ## print tile
 .print_title()

 error.code <- check.inputs.mlam(ncpu=ncpu ,  colname.trait=trait, 
                     map=map, pheno=pheno, geno=geno, Zmat=Zmat, lambda=lambda )
 if(error.code){
   message("\n The Eagle function AM has terminated with errors.\n")
   return(NULL)
 }


 ## checking if map is present. If not, generate a fake map. 
 if(is.null(map)){
   if(!quiet ){
     message(" Map file has not been supplied. An artificial map is being created but this map is not used in the analysis. \n")
     message(" It is only used for the reporting of results. \n")
   }
   ## map has not been supplied. Create own map
   map <- data.frame(SNP=paste("M", 1:geno[["dim_of_M"]][2], sep=""), 
                     Chr=rep(1, geno[["dim_of_M"]][2]), 
                     Pos=1:geno[["dim_of_M"]][2])
  }

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
          message(" fformula has ", length(fformula), " separate terms. It should be of the form  x1 + x2 + x3  \n") 
          message("\n AM has terminated with errors.\n")
          return(NULL)
      }
   } else {
    ## problem: formula should not contain ~
    message(" Only the terms on the right-hand side of the formula should be specified. \n")
    message(" Please remove the ~ from the formula. The formula should be of the form x1 + x2 + x3 \n")
    message("\n AM has terminated with errors.\n")
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
      message("\n  AM has terminated with errors.\n ")
      return(NULL)
   }
  }
 }



 
 ## check for NA's in trait 
 indxNA_pheno <- check.for.NA.in.trait(trait=pheno[[trait]])

 ## set NA's in trait to 0. Extra factor mv will be
 ## fitted by the build_design_matrix function
 if(length(indxNA_pheno) > 0){
       pheno[indxNA_pheno , trait ] <- 0 
 }


 ## build design matrix currentX
 currentX <- .build_design_matrix(pheno=pheno, fformula=fformula, quiet=quiet, 
                 indxNA_pheno=indxNA_pheno)
 currentX <- as.matrix(currentX)

 ## check currentX for solve(crossprod(X, X)) singularity
   chck <- tryCatch({ans <- solve(crossprod(currentX, currentX))},
           error = function(err){
            return(TRUE)
           })

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

  
 ## Initialization
 continue <- TRUE
 itnum <- 1
 outlierstat <- list()
 while(continue){
     message("\n\n Iteration " , itnum, ": Searching for most significant marker-trait association\n\n")
     currentX <- constructX(Zmat=Zmat, fnameMt=geno[["tmpMt"]], currentX=currentX, loci_indx=new_selected_locus,
                          dim_of_Mt=geno[["dim_of_Mt"]],
                          map=map )  
     currentX <- as.matrix(currentX)

     ## calculate Ve and Vg
     Args <- list(geno=geno,
                    ncpu=ncpu,selected_loci=selected_loci,
                    quiet=quiet)

    if(itnum==1){
       if(!quiet)
           message(" quiet=FALSE: calculating M %*% M^t. \n")

       gc()

       MMt <- .calcMMt( geno=geno, ncpu=ncpu, selected_loci= selected_loci, quiet=quiet)
       invMMt <- chol2inv(chol(MMt))   ## doesn't use GPU
       MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=quiet )
       if(is.null(Zmat)){
            eig.L <- emma.eigen.L.wo.Z(MMt )
        } else  {
            eig.L <- emma.eigen.L.w.Z(Zmat, MMt)
        }
        if(!quiet)
            doquiet(dat=MMt, num_markers=5 , lab="M%*%M^t")

     }  ## if itnum==1  


 
     if(!quiet){
        message(" Calculating variance components for multiple-locus model. \n")
     }


     vc <- .calcVC(trait=pheno[, trait ], Zmat=Zmat, currentX=as.matrix(currentX), MMt=MMt, 
                          eig.L=eig.L) 

     gc()
     best_ve <- vc[["ve"]]
     best_vg <- vc[["vg"]]

    
           start <- Sys.time()
     new_extBIC <- .calc_extBIC(vc$ML , pheno[, trait ], currentX, geno, 
                       numberSNPselected=(itnum-1), quiet, lambda) 
           end <- Sys.time()
           #print(c(" .calc_extBIC  ", end-start))

     gc()

     ## set vector extBIC
     extBIC <- c(extBIC, new_extBIC)

     ## Print findings to screen
    .print_results(itnum, selected_loci, map,  extBIC)
  

    ## select new locus if fixit 
    if (fixit){
       if (itnum <= maxit){
           ## find QTL

           ARgs <- list(MMt_sqrt_and_sqrtinv=MMt_sqrt_and_sqrtinv, Zmat=Zmat, geno=geno,availmemGb=geno[["availmemGb"]], selected_loci=selected_loci,
                 MMt=MMt, invMMt=invMMt, best_ve=best_ve, best_vg=best_vg, currentX=as.matrix(currentX),
                 ncpu=ncpu, quiet=quiet, trait=pheno[, trait ], itnum=itnum)


          ## new_selected_locus <- do.call(.find_qtl, ARgs)  ## memory blowing up here !!!! 
          fq <-  do.call(.find_qtl, ARgs)  ## memory blowing up here !!!!

          new_selected_locus <- fq[["orig_indx"]]
          outlierstat[[itnum]] <- fq[["outlierstat"]]

          #print("end")
          gc()
          selected_loci <- c(selected_loci, new_selected_locus)
       } else {
         continue <- FALSE
       }
    } else {

        ## Select new locus if extBIC is still decreasing 
        if(which(extBIC==min(extBIC))==length(extBIC) ){  ## new way of stoppint based on extBIC only

           ## find QTL
           ARgs <- list(MMt_sqrt_and_sqrtinv=MMt_sqrt_and_sqrtinv, Zmat=Zmat, geno=geno,availmemGb=geno[["availmemGb"]], selected_loci=selected_loci,
                 MMt=MMt, invMMt=invMMt, best_ve=best_ve, best_vg=best_vg, currentX=as.matrix(currentX),
                 ncpu=ncpu, quiet=quiet, trait=pheno[, trait ], itnum=itnum)
          start <- Sys.time() 
          fq <-  do.call(.find_qtl, ARgs)  ## memory blowing up here !!!!
          end <- Sys.time() 
           #print(c(" .find_qtl =  ", end-start))


          new_selected_locus <- fq[["orig_indx"]]
          outlierstat[[itnum]] <- fq[["outlierstat"]]
          #print("end")
          gc()
          selected_loci <- c(selected_loci, new_selected_locus)

       }  else {
            ## terminate while loop, 
            continue <- FALSE
       }  ## end if else

   } ## end if fixit
   itnum <- itnum + 1
   ## alternate stopping rule - if maxit has been exceeded.
    if(itnum > maxit){
         continue <- FALSE 
         ## ==> .print_header()
         ## need to remove the last selected locus since we don't go on and calculate its H and extBIC 
         ## under this new model. 
         ##===> .print_final(selected_loci[-length(selected_loci)], map, extBIC, lambda)
         ##===>sigres <- .form_results(traitname=trait, 
         ##                        trait=pheno[, trait ], 
         ##                        selected_loci=selected_loci[-length(selected_loci)], 
         ##                        fformula=fformula,
         ##                        indxNA_pheno=indxNA_pheno,
         ##                        ncpu=ncpu,  availmemGb=geno[["availmemGb"]],
         ##                         quiet=quiet, extBIC=extBIC, lambda=lambda, 
         ##                        geno=geno, pheno=pheno, map=map, Zmat=Zmat, outlierstat=outlierstat) 
    }
 
  }  ## end while continue

if( itnum > maxit){
   .print_header()
   .print_final(selected_loci[-length(selected_loci)], map,  extBIC, lambda)
message("The maximum number of iterations ", maxit, " has been reached. ")
message("To increase the maximum number of iterations, adjust ")
message(" the maxit parameter in the AM function call. \n ")
    sigres <- .form_results(traitname=trait, trait=pheno[, trait ], selected_loci[-length(selected_loci)],   fformula, 
                     indxNA_pheno,  ncpu, geno[["availmemGb"]], quiet,  extBIC, lambda,
                     geno, pheno, map, Zmat, outlierstat )   

} else {
    ## remove last selected_loci as for this locus, the extBIC went up
    if(length(selected_loci)>1){
        .print_header()
        .print_final(selected_loci[-length(selected_loci)], 
                     map, 
                     extBIC[-length(selected_loci)], lambda )
        sigres <- .form_results(traitname=trait, pheno[, trait ], selected_loci[-length(selected_loci)],   fformula, 
                         indxNA_pheno,  ncpu,  geno[["availmemGb"]], quiet, 
                         extBIC[-length(selected_loci)], lambda, 
                         geno, pheno, map, Zmat, outlierstat )   
    } else {
        .print_header()
        .print_final(selected_loci, map, extBIC, lambda )
        sigres <- .form_results(traitname=trait, pheno[, trait ], selected_loci,   fformula, 
                         indxNA_pheno,  ncpu,  geno[["availmemGb"]], quiet, extBIC, lambda, 
                         geno, pheno, map, Zmat, outlierstat )   
   }  ## end inner  if(length(selected_locus)>1)
}  ## end if( itnum > maxit)

return( sigres )

} ## end AM


computer <- new.env()
computer$ngpu <- 0






