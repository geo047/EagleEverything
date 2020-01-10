create.ascii  <- function(file_genotype=NULL,  type="text", AA=NULL, AB=NULL, BB=NULL,
                         availmemGb=8, dim_of_M=NULL, quiet=TRUE, missing=NULL){
 ## an Rcpp function to create the no-space file of the genotype data M and Mt
 ## from marker data. The marker data may be from an ASCII file or PLINK ped file.
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## AA, AB, BB       numeric codes for associated genotypes in marker genotype file
 ## availmemGb     available memory for conversion to reformatted file
 ## dim_of_M             row, column dimensions of M.  
 ## type            where file type is text or PLINK

 if(.Platform$OS.type == "unix") {
       ##tmpMfile <- paste(dirname(file_genotype), "/", "M.tmp", sep="")
       ##tmpMtfile <- paste(dirname(file_genotype), "/", "Mt.tmp", sep="")
       tmpMfile <- paste(tempdir() , "/", "M.tmp", sep="")
       tmpMtfile <- paste(tempdir()  , "/", "Mt.tmp", sep="")
 } else {
       ##tmpMfile <- paste(dirname(file_genotype), "\\", "M.tmp", sep="")
       ##tmpMtfile <- paste(dirname(file_genotype), "\\", "Mt.tmp", sep="")
       tmpMfile <- paste(tempdir() , "\\", "M.tmp", sep="")
       tmpMtfile <- paste(tempdir() , "\\", "Mt.tmp", sep="")
 }



if (type=="text"){
    ## text genotype file
    if(!is.null(missing)) {
        missing <- as.character(missing)
    } else {
      missing <- "NA"
    }
    it_worked <- createM_ASCII_rcpp(f_name = file_genotype, type=type ,  f_name_tmp = tmpMfile, AA = AA, AB = AB, BB = BB,
               max_memory_in_Gbytes=availmemGb,  dims = dim_of_M ,
               quiet = quiet, missing=missing, message=message)
    if(!it_worked) #  creation of ASCII file has failed 
       return(FALSE)

    message(" \n Taking transpose of marker data and writing untransposed and transposed data to disc ... \n") 
    createMt_ASCII_rcpp(f_name = tmpMfile, f_name_tmp = tmpMtfile,   type=type,
                  max_memory_in_Gbytes=availmemGb,  dims = dim_of_M, quiet = quiet , message=message )
    message("\n  Writing of marker data to disc is complete ... \n")
} else {
    ## PLINK ped file
    ## using -9 to indicate missing/null genotypes
    ncol  <- dim_of_M[2]
    dim_of_M[2] <- 2*dim_of_M[2] + 6  ## number of cols in a PLINK file
    it_worked <- createM_ASCII_rcpp(f_name = file_genotype, type=type,  f_name_tmp = tmpMfile, AA ="-9", AB = "-9", BB = "-9",
               max_memory_in_Gbytes=availmemGb,  dims = dim_of_M , quiet = quiet,
                missing="NA", message=message)
     if(!it_worked) #  creation of ASCII file has failed 
       return(FALSE)


    dim_of_M[2] <- ncol ## setting back to number of cols in no-space ASCII file
    message(" \n Taking transpose of marker data and writing untransposed and transposed data to disc ... \n") 
    createMt_ASCII_rcpp(f_name = tmpMfile, f_name_tmp = tmpMtfile,    type=type,
                  max_memory_in_Gbytes=availmemGb,  dims = dim_of_M, quiet = quiet , message=message )
    message(" \n Writing of marker data to disc is complete ... \n")

}  ## end if else type

 return(TRUE)

}



