#' @title Calculate false positive rate (FPR) for AM
#' @description A permutation approach for calculating the false positive rate for \code{\link{AM}} for a given \code{gamma} value. 
#' 
#' @param gamma     a value between 0 and 1. Values closer to 1 decrease false positive rate.  This parameter must be specified.  
#' @param  numreps  the number of replicates upon which to base the calculation of the false 
#'                   positive rate. We have found 100 replicates to be sufficient.  
#' @param trait  the name of the column in the phenotype data file that contains the trait data. The name is case sensitive and must match exactly the column name in the phenotype data file.  This parameter must be specified. 
#' @param fformula   the right hand side formula for the fixed effects part of the model. 
#' @param availmemGb a numeric value. It specifies the amount of available memory (in Gigabytes). 
#' This should be set to the maximum practical value of available memory for the analysis. 
#' @param geno   the R  object obtained from running \code{\link{ReadMarker}}. This must be specified. 
#' @param pheno  the R  object  obtained  from running \code{\link{ReadPheno}}. This must be specified.
#' @param map   the R object obtained from running \code{\link{ReadMap}}. If not specified, a generic map will 
#'              be assumed. 
#' @param Zmat     the R object obtained from running \code{\link{ReadZmat}}. If not specified, an identity matrix will be assumed. 
#' @param ncpu a integer  value for the number of CPU that are available for distributed computing.  The default is to determine the number of CPU automatically. 
#' @param ngpu   a integer value for the number of gpu available for computation.  The default
#'               is to assume there are no gpu available.  
#'               This option has not yet been implemented.
#' @param seed  a integer value for the starting seed for the permutations. 
#' 
#' @details
#'
#' This function is useful for  calculating  the false positive rate  of \code{\link{AM}} for a given value of \code{gamma}. We do this empirically, via permutation. 
#' This function is also of use when interested in finding the  value of \code{gamma}  corresponding to a certain false positive rate. Ideally, this would be done automatically 
#' but currently, it is up to the user to implement their own optimisation strategy. It is on our 'to do list' to revise this function so that 
#'  a binary search is performed to find the 'best' \code{gamma} value. 
#'
#'
#'
#' Permutation is performed as follows. First, the null model is fitted to the trait data. 
#' This is the model that contains all the fixed effects and an error. 
#' Second,  residuals are  obtained 
#' from the fitted null model. Third, the residuals are permuted \code{numreps} times. Fourth, for each permutation, 
#' we treat the residuals as the trait and perform multiple-locus association mapping, conditional on the \code{gamma}
#' value that has been specified. Any findings are false positives. Fifth, once the \code{numreps} analyses of the permuted 
#' residuals have been performed, we sum the total number of false positives across the replicate analyses, and 
#' divide by \code{numreps}.  We now have an estimate of the false positive rate (for the specified trait, 
#' fixed effects model \code{fformula} and \code{gamma} value). 
#'
#'
#' @seealso \code{\link{AM}}
#' @return
#' A list with the following components:
#' \describe{
#'\item{numreps}{the number of permutations performed.}
#'\item{gamma}{gamma value that was used to calculate the false positive rate}
#'\item{falsepos}{the false positive rate, given the gamma value, calculated from 
#' analyses of the 'numreps'  permutations.}
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
#'  # Find the false positive rate for AM when gamma is set to the
#'  # value 0.9
#'  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'  
#'   falseposrate <- FPR4AM(trait = 'y',
#'                 gamma= 0.9,
#'                 fformula=c('cov1+cov2'),
#'                 map = map_obj,
#'                 pheno = pheno_obj,
#'                 geno = geno_obj) 
#'
#' }
#'
#'
FPR4AM <- function(
               trait=trait,
               gamma = NULL,
               numreps = 100,
               fformula  = NULL,
               availmemGb=8,
               geno=NULL,
               pheno=NULL,
               map = NULL,
               Zmat = NULL,
               ncpu=detectCores(),
               ngpu=0,
               seed=101
               ){

  set.seed(seed)
 # need some checks in here ... 
error.code <- check.inputs.mlam(ncpu=ncpu , availmemGb=availmemGb, colname.trait=trait,
                     map=map, pheno=pheno, geno=geno, Zmat=Zmat )
 if(error.code){
   message("\n The Eagle function CalculateFDR has terminated with errors.\n")
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

 # fit null model and find residuals
 pheno[, "trait"] <- pheno[, trait]  ## name change
 if(is.null(fformula)){
   ff <- as.formula("trait ~ 1")
 } else {
    # turn the formula back into a character
    fformula <- as.character(fformula)[2]  # gets rid of the ~
    ff <- paste("trait ~ ", fformula, sep=" ")
    ff <- as.formula(ff)
 }

 mod <- lm( ff , data=pheno)
 res <- residuals(mod)
 # if NA's in trait or fixed effects, then residuals are not calculated. 
 # Need to account for this. 
 indx <- as.numeric(names(res))
 pheno[, "residuals"] <- NA  # initialize
 pheno[indx, "residuals"] <- res

am_res <- list()
for (ii in 1:numreps){
  ## permute residuals and perform analysis (only permute non-NA values)
  indx <- which(!is.na(pheno[, "residuals"]))
  sperm <- pheno[indx, "residuals"]
  sperm <- sperm[sample(1:length(sperm), length(sperm), FALSE)]
  pheno[ , "residuals"] <- NA
  pheno[indx, "residuals"] <- sperm
  #pheno[, "residuals"] <- pheno[sample(1:length(res), length(res), FALSE), "residuals" ]
  
   ## Perform AM+ analysis
   argu <- list(
               trait="residuals",
               gamma = gamma, 
               availmemGb=availmemGb,
               geno=geno,
               pheno=pheno,
               map = map,
               ncpu=ncpu, 
               ngpu=ngpu)

   am_res[[ii]] <- do.call( AM, argu)

   ## running total
   numres <- rep(0, ii)
   for(jj in 1:ii){
       indx <- am_res[[jj]]$Mrk
       indx <- indx[!is.na(indx)]
       if(length(indx) > 0)
           numres[jj] <- length(indx)
       }

       falsepos <- sum(numres)/ii

       cat(" Iteration ", ii, "\n")
       #cat(" Gamma = ", gamma, "\n")
       #cat(" False positive (FP) rate is ", round( falsepos, 3)   , " with ", sum(numres), "FPs already found \n")


  
}  ## end for ii in 1:numreps

numres <- rep(0, numreps)
for(ii in 1:numreps){
  indx <- am_res[[ii]]$Mrk
  indx <- indx[!is.na(indx)]
  if(length(indx) > 0)
     numres[ii] <- length(indx)
}

  falsepos <- sum(numres)/numreps

  cat(" For a gamma of ", gamma, " the false positive rate is ", 
           round( falsepos, 3)   , " \n")
  cat(" If the false positive rate is too high (low), increase (decrease) the value of gamma. \n")


 return(list(numreps=numreps, gamma=gamma, falsepos=falsepos    ))


}





