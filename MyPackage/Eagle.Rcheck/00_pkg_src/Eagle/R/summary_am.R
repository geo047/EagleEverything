#' @title Summary of multiple locus association mapping results
#' @description    A summary function that provides additional information on the significant 
#'     marker-trait associations found by \code{\link{AM}}
#' @param  AMobj  the (list) object obtained from running \code{\link{AM}}. 
#' @details
#'
#' \code{SummaryAM} produces two  tables, an overall summary table and a table  of results with 
#' the  p-value for each 
#' fixed effect in the final model.  
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
#'                            geno = geno_obj)
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

 ngpu <- AMobj$ngpu


 if(is.null(AMobj)){
    message(" SummaryAM function requires AMobj object to be specified. This object is obtained by running AM().")
    return(NULL)
    }
 if(!is.list(AMobj)){
    message(" SummaryAM function requires AMobj object to be a list object.")
    return(NULL)
   }




## Table 1
# create list with necessary components 
lst <- list( "ncpu" = AMobj[["ncpu"]],
             "memory" = AMobj[["availmemGb"]],
             "numsnp"=nrow(AMobj[["map"]]),
            "numsamples"=nrow(AMobj[["pheno"]]),
            "traitname"=AMobj[["traitname"]], 
            "fformula" = AMobj[["fformula"]],
            "nummissingpheno" = AMobj[["indxNA_pheno"]],
            "numSigSNP"=AMobj[["Mrk"]],
            "gamma" = AMobj[["gamma"]] 
         )

if (is.null(lst[["fformula"]])){
 lst[["fformula"]] <- "intercept only"
} else {
  lst[["fformula"]] <- as.character(AMobj[["fformula"]])[2]
} 

if (is.null(lst[["nummissingpheno"]])){
lst[["nummissingpheno"]] <- 0
} else {
lst[["nummissingpheno"]] <- length(lst[["nummissingpheno"]])
}
if (length(lst[["numSigSNP"]]) == 1) {
  lst[["numSigSNP"]] <- 0
} else {
  lst[["numSigSNP"]] <- length(lst[["numSigSNP"]]) - 1
}

message("\n\n Table 1: Summary Information \n   ")
message(  sprintf("%50s", "--------------------------------------------------------" ))
message( sprintf("%-40s  %-10s", "Number cpu: ", lst[["ncpu"]] ))
message( sprintf("%-40s  %-10s", "Max memory (Gb): ", lst[["memory"]] ))
message( sprintf("%-40s  %-10s", "Number of samples: ", lst[["numsamples"]] ))
message( sprintf("%-40s  %-10s", "Number of snp: ", lst[["numsnp"]] ))
message( sprintf("%-40s  %-10s", "Trait name: ", lst[["traitname"]] ))
message( sprintf("%-40s  %-30s", "Fixed model: ", lst[["fformula"]] ))
message( sprintf("%-40s  %-10s", "Number samples missing obs:", lst[["nummissingpheno"]] ))
message( sprintf("%-40s  %-10s", "Number significant snp-trait assocs:", lst[["numSigSNP"]] ))
message( sprintf("%-40s  %-4s", "Gamma value for extBIC: ", round(lst[["gamma"]],2)  ))
message(  sprintf("%50s", "--------------------------------------------------------" ))
message("\n\n")

# create data frame of summary information for use in shiny app
infodf <- data.frame("description"= c("Number cpu", "Max memory (Gb)", "Number of samples", "Number of snp", 
                                    "Trait name", "Fixed model", "Number samples missing obs", 
                                    "Number significant snp-trait assocs", "Gamma value for extBIC"  )   ,
                   "value" = c(lst[["ncpu"]], lst[["memory"]], lst[["numsamples"]], 
                     lst[["numsnp"]], lst[["traitname"]], lst[["fformula"]], lst[["nummissingpheno"]], 
                    lst[["numSigSNP"]], round(lst[["gamma"]],2)   )
         )


## Table II
 ## check to make sure that null model is not being supplied
 if (length(AMobj$Mrk)==1){
   message(" No significant marker-trait associations have been found by AM. \n")
   message(" No p-values to report \n")
   return(invisible(lst) )
 }

  ## build environmental effects design matrix
  baseX <- .build_design_matrix(pheno=AMobj$pheno,  
                                    fformula=AMobj$fformula,
                                   quiet=AMobj$quiet, indxNA_pheno=AMobj$indxNA_pheno)

  ## add genetic marker effects 
  fullX <- baseX
  for(loc in AMobj$Indx){
           
           fullX <- constructX(Zmat=AMobj$Zmat, fnameMt=AMobj$geno[["asciifileMt"]], 
                              currentX=fullX, loci_indx=loc,
                               dim_of_Mt=AMobj$geno[["dim_of_ascii_Mt"]],
                                map=AMobj$map)
  }  ## end for loc

  ## calculate MMt
  MMt <- .calcMMt(AMobj$geno,  AMobj$ncpu, AMobj$Indx, AMobj$quiet)

  ## calculate variance components of LMM

   Args <- list(y= AMobj$trait  , X= fullX , Z=AMobj$Zmat, K=MMt, ngpu=ngpu)

  #eR <- emma.REMLE(y=AMobj$trait, X= fullX , K=MMt, llim=-100,ulim=100)
  eR <-  do.call(emma.MLE, Args)






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

message("\n\n Table 2: Significance of Effects in Final Model \n   ")

  message(sprintf("%20s %6s %15s      %8s", "",   "Df", "Wald statstic" , "Pr(Chisq)"))
  for(ii in varnames )
  {
      indx <- which(varnames == ii)
      message(sprintf("%20s %6i     %10.2f       %.3E",
         ii,   df[indx], W[indx], pval[indx ]))
  }  ## end for ii
 message("\n\n\n")
df_size <- data.frame("Effects"=varnames,  "Df"=as.character(df),   
                      "Wald statistic"=as.character(round(W,2)),        
                      "Pr(Chisq)"=pval, check.names=FALSE)


 ##----------------------------------------------------------------------- 
 ## Variance explained - based on Sun et al. (2010). Heredity 105:333-340
 ##----------------------------------------------------------------------- 
 # 13 Nov, 2019
 # Put on hold - not confident in the results. 

# MMt <- MMt/max(MMt) + 0.05 * diag(nrow(MMt))  
# # base model
# # Note - baseX = fullX so okay to leave Args as is
# #basemod <- emma.MLE(y=AMobj$trait, X=baseX, K=MMt, llim=-100,ulim=100)
# Args[["X"]] <- baseX
#
# basemod  <-  do.call(emma.MLE, Args)
#
#
#
# base_logML <- basemod$ML
# # full model
#  df_R <- NULL
#  fullX <- baseX
#  message(" Table 2: Proportion of phenotypic variance explained by the ")
#  message("          model. Marker loci, which were found by AM(), are added")
#  message("          a SNP at a time.\n ")
#  message(sprintf("    %18s          %10s ", "SNP", "Proportion"))
#  for(loc in AMobj$Indx[-1]){
#               fullX <- constructX(Zmat=AMobj$Zmat, fnameMt=AMobj$geno[["asciifileMt"]],
#                                currentX=fullX, loci_indx=loc,
#                               dim_of_Mt=AMobj$geno[["dim_of_ascii_Mt"]],
#                               map=AMobj$map)
#        # fullmod <- emma.MLE(y=AMobj$trait, X=fullX, K=MMt, llim=-100,ulim=100)
#        ## calculate variance components of LMM
#        Args[["X"]]  <-  fullX 
#        fullmod  <-  do.call(emma.MLE, Args)
#        full_logML <- fullmod$ML
#        print(fullX[1:5,])
#        cat(" Base logML = ", base_logML , "\n")
#        cat(" Full logML = ", full_logML, "\n")
#        Rsq <- 1 - exp(-2/length(AMobj$trait) * (full_logML - base_logML))
#        print(Rsq)
#        message(sprintf("  %+20s               %.3f",  paste("+ ",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])), Rsq))
#        df_R <- rbind.data.frame(df_R, data.frame("Marker_name"=paste("+",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])),
#                                                  "Prop_var_explained"=Rsq))
#   }  ## end for loc
#  colnames(df_R) <- c("SNP", "Prop var explained")
#
  res <- list()
  res[["pvalue"]] <- pval
  res[["size"]] <- df_size
  res[["Waldstat"]] <- W
  res[["df"]] <- df
  res[["summarylist"]] <- infodf

  
  return(invisible(res))
}
