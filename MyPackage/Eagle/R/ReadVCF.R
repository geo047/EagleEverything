
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
#'  VCF is a tab separated text file containing  meta-information lines, a header line, and then data 
#' lines. The data lines  contain information about a position in the genome.
#'
#' It is assumed that genotype information has been recorded on samples for each position. 
#'
#' Loci with more than two alleles will be removed automatically. 
#' 
#' Eagle will only accept a single (uncompressed) vcf file. If chromosomal information has been recorded in separate 
#' vcf files, these files need to be merged into a single vcf file.  This can be done by using the 
#'  BCFtools utility set with command line "bcftools concat".
#'
#'
#' @return  To allow Eagle to handle data larger than the memory capacity of a machine, \code{ReadVCF} doesn't load 
#' the marker data into memory. Instead, it 
#' writes a reformatted version of the marker data, and its transpose, to the harddrive. These two files
#' are only temporary, being removed at the end of the R session. 
#' The object returned by
#' \code{ReadVCF} is a list object with the elements \code{tmpM} , \code{tmpMt}, and \code{dim_of_M}  
#' which is the full file name (name and path)  of the reformatted file for the marker  data,  the full file name of the reformatted file 
#' for the transpose of the marker  data,  and a 2 element vector with the first element the number of individuals and the second 
#' element the number of marker loci. 
#' 
#'
#' @examples
#'   #
#'   # Read in the genotype data contained in the vcf file geno.vcf
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file('extdata', 'geno.vcf', package='Eagle')
#'   # 
#'   # The full path and name of the file is
#'   print(complete.name)
#'   
#'   # The file contains 5 marker loci recorded on 3 individuals
#'   # Two of the loci contain multiple alleles and are removed. 
#'   # A summary of the file is printed once the file has been read.
#'   geno_obj <- ReadVCF(filename=complete.name, availmemGb=4) 
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
               "dim_of_M" = liststr[["dim_of_M"]],
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



