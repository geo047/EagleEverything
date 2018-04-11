#' @title Calculate false positive rate for extBIC 
#' @description Used for fine tuning the gamma value in the extBIC.
#' 
#' @param gamma     a value between 0 and 1 for the extBIC. 
#' Values closer to 0 makes Eagle more prone to finding a false positive but increases the 
#' chance of finding a true positive. Values closer to 1 makes Eagle more conservative, with 
#' less chance of finding a false positive  but also less chance of finding true 
#' positives. 
#' @param  numreps  the number of replicates upon which to base the calculation of the false 
#'                   positive rate. 
#' @param trait  the name of the column in the phenotype data file that contains the trait data. The name is case sensitive and must match exactly the column name in the phenotype data file. 
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
#' 
#' @details
#'
#' This function performs  \code{numreps} analyses of permuted data, from which the 
#' false positive rate for a given \code{gamma} value can be calculated, empirically. 
#' This function is of use for setting gamma to its optimal value. 
#'
#' 
#' Background: Eagle uses the extended BIC (extBIC) in which to decide if it should keep looking for more 
#' significant associations or to stop.  The conservativeness of  extBIC is adjusted via the 
#' gamma parameter. Values closer to 0 make the extBIC less conservative. Values closer to 1 
#' make the extBIC more conservative. 
#' If a gamma value is not specified in \code{\link{AM}}, then a default value is assigned 
#' based on sample size. This default value is assigned according to a heuristic procedure. 
#'
#'
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
#'
CalculateFPR <- function(
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

  for(ii in 1:numreps){
    pheno[, "trait"] <- pheno[, trait]  ## name change
    ff <- paste("trait ~ ", fformula, sep=" ")
    ff <- as.formula(ff)
    mod <- lm( ff , data=pheno)
    res <- as.vector(residuals(mod))
    pheno[, "residuals"] <- res[sample(1:length(res), length(res), FALSE)]
  }


am_res <- list()
for (ii in 1:numreps){
  ## find the residuals
  pheno[, "trait"] <- pheno[, trait]  ## name change
  ff <- paste("trait ~ ", fformula, sep=" ")
  ff <- as.formula(ff)
  mod <- lm( ff , data=pheno)
  res <- as.vector(residuals(mod))
  pheno[, "residuals"] <- res[sample(1:length(res), length(res), FALSE)]
  

 
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
  cat(" Gamma = ", gamma, "\n")
  cat(" False positive rate is ", round( falsepos, 3)   , " with ", sum(numres), "FPs already found \n")


  
}

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



 return(list(numreps=numreps, gamma=gamma, falsepos=falsepos    ))


}





