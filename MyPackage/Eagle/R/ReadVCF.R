
#' @title Read  vcf file
#' 
#' @description
#' A function for reading in marker data in vcf format. 
#' @param filename contains the name of the vcf  file. The file name needs to be in quotes. If the file is not in the working directory, then the full path 
#' to the file is required.
#' @param availmemGb a numeric value. It specifies the amount of available memory (in Gigabytes). 
#'         This should be set to be as large as possible for best performance.   
#' @param  quiet      a logical value. If set to \code{TRUE}, additional runtime output is printed. 
#' @details
#' 
#'
#' @return  To allow Eagle to handle data larger than the memory capacity of a machine, \code{ReadVCF} doesn't load 
#' the marker data into memory. Instead, it creates a reformatted file of the marker data and its transpose. The object returned by
#' \code{ReadVCF} is a list object with the elements \code{tmpM} , \code{tmpMt}, and \code{dim_of_M}  
#' which is the full file name (name and path)  of the reformatted file for the marker  data,  the full file name of the reformatted file 
#' for the transpose of the marker  data,  and a 2 element vector with the first element the number of individuals and the second 
#' element the number of marker loci. 
#' 
#'
#' @examples
#'   #
#'   # Read in the genotype data contained in the text file geno.txt
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file('extdata', 'geno.txt', package='Eagle')
#'   # 
#'   # The full path and name of the file is
#'   print(complete.name)
#'   
#'   # Here, 0 values are being treated as genotype AA,
#'   # 1 values are being treated as genotype AB, 
#'   # and 2 values are being treated as genotype BB. 
#'   # 4 gigabytes of memory has been specified. 
#'   # The file is space separated with the rows the individuals
#'   # and the columns the snp loci.
#'   geno_obj <- ReadVCF(filename=complete.name, type='text', AA=0, AB=1, BB=2, availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
#'   #--------------------------------
#'   #  Example 2
#'   #-------------------------------
#'   #
#'   # Read in the allelic data contained in the PLINK ped file geno.ped
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file('extdata', 'geno.ped', package='Eagle')
#'
#'   # 
#'   # The full path and name of the file is
#'   print(complete.name)
#'   
#'   # Here,  the first 6 columns are being ignored and the allelic 
#'   # information in columns 7 -  10002 is being converted into a reformatted file. 
#'   # 4 gigabytes of memory has been specified. 
#'   # The file is space separated with the rows the individuals
#'   # and the columns the snp loci.
#'   geno_obj <- ReadVCF(filename=complete.name, type='PLINK', availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
#'
ReadVCF <- function( filename=NULL, availmemGb=16, quiet=TRUE ){


 if (nargs() == 0){
    ## checking that function has arguments
    message(" Please supply arguments to function \n")
    return(NULL)
 }

       
 ## checking if a VCF file has been specified. 
 if (is.null(filename)){
            message(" The name of the vcf file is missing.")
            message(" ReadVCF has terminated with errors.")
            return(NULL)
 }
 if (!file.exists(fullpath(filename) )){
            message(" The vcf file ", filename, " could not be found. ")
            message(" ReadVCF has terminated with errors ")
            return(NULL)
 }

 ## Rcpp function to create binary packed M and Mt file 
 #dim_of_M <- create.vcf.bin(file_genotype=fullpath(filename), availmemGb=availmemGb,  quiet=quiet  )
 liststr <- create.vcf.bin(file_genotype=fullpath(filename), availmemGb=availmemGb,  quiet=quiet  )


    if(.Platform$OS.type == "unix") {
       binfileM <- paste(tempdir(), "/", "M.bin", sep="")
       binfileMt <- paste(tempdir(), "/", "Mt.bin", sep="")
     } else {
       binfileM <- paste(tempdir()  , "\\", "M.bin", sep="")
       binfileMt <- paste(tempdir() , "\\", "Mt.bin", sep="")
     }


 geno <- list("tmpM"=binfileM, "tmpMt"=binfileMt,
               "dim_of__M" = liststr[["dim_of_M"]],
               "dim_of_Mt" = c( liststr[["dim_of_M"]][2], liststr[["dim_of_M"]][1]),
               "availmemGb" = availmemGb, 
               "map" = liststr[["map"]] )



  if(.Platform$OS.type == "unix") {
       RDatafile <- paste(tempdir() , "/", "M.RData", sep="")
  } else {
       RDatafile <- paste( tempdir() , "\\", "M.RData", sep="")
  }




  save(geno, file=RDatafile)
  ## create M.Rdata file in current directory
  return(geno)


}  ## end function call ReadVCF



