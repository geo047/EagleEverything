#' @title Summary of multiple locus association mapping results
#' @description    A summary function that provides additional information on the significant 
#'     marker-trait associations found by \code{\link{AM}}
#' @param  AMobj  the (list) object obtained from running \code{\link{AM}}. 
#' @details
#'
#' \code{SummaryAM} produces two tables of results. First, a table of results is produced with 
#' the additive effect size and p-value for each 
#' fixed effect in the final model.  Second, a table of results is produced with the 
#' proportion of phenotypes variance explained by  the different multiple-locus models. Each row 
#' in this table is the proportion of phenotypic variance explained (Sun et al. 2010) after the marker locus has been added to the 
#' multiple locus model. 
#' @references  Sun G., Zhu C., Kramer  MH., Yang S-S., et al. 2010. Variation explained in mixed model association 
#' mapping. Heredity 105, 330-340. 
#' @examples
#'  \dontrun{
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #---------------
#'   # read the map 
#'   #---------------
#'   #
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   #------------------
#'   # read marker data
#'   #------------------
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8 gigabytes
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   #----------------------
#'   # read phenotype data
#'   #-----------------------
#'
#'   # Read in a plain text file with data on a single trait and two fixed effects
#'   # The first row of the text file contains the column names y, cov1, and cov2. 
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   #-------------------------------------------------------
#'   # Perform multiple-locus genome-wide association mapping 
#'   #-------------------------------------------------------                   
#'   res <- AM(trait = 'y',
#'                            fformula=c("cov1 + cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
#'
#'   #-----------------------------------------
#'   # Produce additional summary information 
#'   #------------------------------------------
#'
#'   SummaryAM(AMobj=res)
#'  }
#'
#' 
#' 
#' @seealso \code{\link{AM}}
#'
SummaryAM <- function(AMobj=NULL)
{

 if(is.null(AMobj)){
    message(" SummaryAM function requires AMobj object to be specified.")
    return(NULL)
    }
 if(!is.list(AMobj)){
    message(" SummaryAM function requires AMobj object to be a list object.")
    return(NULL)
   }

 ## check to make sure that null model is not being supplied
 if (length(AMobj$Mrk)==1){
   message(" No significant marker-trait associations have been found by AM. \n")
   message(" Nothing to summarize. \n")
   return()
 }

  ## build environmental effects design matrix
  baseX <- .build_design_matrix(pheno=AMobj$pheno,  indxNA=AMobj$indxNA_pheno, 
                                    fformula=AMobj$fformula,
                                   quiet=AMobj$quiet)

  ## add genetic marker effects 
  fullX <- baseX
  for(loc in AMobj$Indx){
           fullX <- constructX(fnameM=AMobj$geno[["asciifileM"]], 
                              currentX=fullX, loci_indx=loc,
                               dim_of_ascii_M=AMobj$geno[["dim_of_ascii_M"]],
                                map=AMobj$map)
  }  ## end for loc

  ## calculate MMt
  MMt <- .calcMMt(AMobj$geno, AMobj$availmemGb, AMobj$ncpu, AMobj$Indx, AMobj$quiet)

  ## calculate variance components of LMM
  eR <- emma.REMLE(y=AMobj$trait, X= fullX , K=MMt, llim=-100,ulim=100)

 ## calculating p values of fixed marker effect via Wald statistic
 mrks <- AMobj$Mrk[-1]  ## its -1 to remove the NA for the null model 
 pval <- vector("numeric", length(colnames(fullX)) )


 H <- calculateH(MMt=MMt, varE=eR$ve, varG=eR$vg, Zmat=AMobj$Zmat, message=message )

# H <-  eR$vg * MMt + eR$ve * diag(1, nrow(MMt))
 Hinv <- try(solve(H))
 beta <- try(solve( t(fullX) %*% Hinv %*% fullX) %*% t(fullX) %*% Hinv %*% matrix(data=AMobj$trait ,ncol=1)   )

# code to work out degrees of freedom and marker names for effects
# in fformula. 
if(is.null(AMobj$fformula)){
  # no fixed effects beside default intercept in model
  df <- 1
  varnames <- "(Intercept)"
} else {
   df_intercept  <- 1  ## intercept only model
   ## get variable names but doesn't include intercept
   nms <-  names(get_all_vars(formula=as.formula(AMobj$fformula), data=AMobj$pheno))

   ## dealing with variables in fformula
   if(length(nms) > 0){
      df_effects <- c(rep(NA, length(nms)))  # intercept + other effects
      # calculate degrees of freedom for variables
      for(ii in 1:length(nms))
      {
  
          if(is.null(levels(AMobj$pheno[, nms[ii]]))){
              df_effects[ii] <- 1   ## its a covariate
          } else {
            df_effects[ii] <- length(levels(AMobj$pheno[, nms[ii]])) - 1  ## for a factor 
          }
      } 
    } ## end if

   df <- df_intercept
   if(length(nms) > 0)
     df <- c(df_intercept, df_effects)

   varnames <- "(Intercept)"
   if(length(nms) > 0)
        varnames <- c(varnames, nms)
} ## end if else

## add markers  to varnames and degrees of freedom

varnames <- c(varnames,  as.character(AMobj$Mrk[!is.na(AMobj$Mrk) ]))
df <- c(df, rep(1, length( AMobj$Mrk[ !is.na(AMobj$Mrk) ] )))


 
pval <- rep(NA, length(varnames))
W <- rep(NA, length(varnames))

for(ii in varnames){
   indx <- grep(ii, colnames(fullX))
   L <- matrix(data=rep(0, length(indx)*ncol(fullX)), byrow=TRUE, nrow=length(indx) )  # r x p matrix
   LL <- diag(length(indx))
   L[,indx] <- LL



  W[which(ii==varnames)]  <- t(L %*% beta) %*%
            solve( L %*% solve(t(fullX) %*% Hinv %*% fullX) %*% t(L) ) %*%
            (L %*% beta)
 pval[which(ii==varnames)] <- 1 - pchisq(W[which(ii==varnames)], length(indx))  ## its not -1 here because fullX already has 1
                                                      ## less factor level
 }

message("\n\n Table 1: Size and Significance of Effects in Final Model \n   ")

  message(sprintf("%20s %10s %6s %15s      %8s", "",  "Size", "Df", "Wald statstic" , "Pr(Chisq)"))
  for(ii in varnames )
  {
      indx <- which(varnames == ii)
      message(sprintf("%20s %10.2f %6i     %10.2f       %.3E",
         ii,  beta[indx], df[indx], W[indx], pval[indx ]))
  }  ## end for ii
 message("\n\n\n")
df_size <- data.frame("Effects"=varnames, "Size"=beta, "Df"=as.character(df),   
                      "Wald statistic"=as.character(round(W,2)),        
                      "Pr(Chisq)"=pval, check.names=FALSE)



 ##----------------------------------------------------------------------- 
 ## Variance explained - based on Sun et al. (2010). Heredity 105:333-340
 ##----------------------------------------------------------------------- 

 MMt <- MMt/max(MMt) + 0.05 * diag(nrow(MMt))  
 # base model
 basemod <- emma.MLE(y=AMobj$trait, X=baseX, K=MMt, llim=-100,ulim=100)
 base_logML <- basemod$ML

 # full model
  df_R <- NULL
  fullX <- baseX
  message(" Table 2: Proportion of phenotypic variance explained by the ")
  message("          model. Marker loci, which were found by AM(), are added")
  message("          a SNP at a time.\n ")
  message(sprintf("    %18s          %10s ", "SNP", "Proportion"))
  for(loc in AMobj$Indx[-1]){
        fullX <- constructX(fnameM=AMobj$geno[["asciifileM"]],
                                currentX=fullX, loci_indx=loc,
                               dim_of_ascii_M=AMobj$geno[["dim_of_ascii_M"]],
                               map=AMobj$map)
        fullmod <- emma.MLE(y=AMobj$trait, X=fullX, K=MMt, llim=-100,ulim=100)
        full_logML <- fullmod$ML
        Rsq <- 1 - exp(-2/nrow(MMt) * (full_logML - base_logML))
        message(sprintf("  %+20s               %.3f",  paste("+ ",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])), Rsq))
        df_R <- rbind.data.frame(df_R, data.frame("Marker_name"=paste("+",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])),
                                                  "Prop_var_explained"=Rsq))
   }  ## end for loc
  colnames(df_R) <- c("SNP", "Prop var explained")

  res <- list()
  res[["pvalue"]] <- pval
  res[["size"]] <- df_size
  res[["Waldstat"]] <- W
  res[["df"]] <- df
  res[["R"]] <- df_R

  
  invisible(res)
}
