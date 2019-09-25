magmaSolve <- function(Xmat , ngpu=1, printInfo=FALSE){


## R frontend to the magma code magma_solve.
## Xmat is assumed to be symmetric.  
## Due to the 32 bit interface problem between R and Magma, the Inv matrix solution from 
## magma_solve is passed back into R via a temporary binary file using the readBin() function. 

 if (nrow(Xmat) != ncol(Xmat)){
   message("\n magmaSolve function needs a square matrix. \n")
   return(0)
 }

 # A temporary file name is assigned for the cpp function to use
 binInvfile <- paste(tempdir() , "/", "Inv.bin", sep="")


 success <- magma_solve(X=Xmat ,  numgpus=ngpu, printInfo=printInfo, fname=binInvfile )
 if(success !=  0)
  {
   message("\n Error: magmaSolve function has failed. A zero value has been returned.  \n")
   return(0)
  } else {

 MaxIntVal <- .Machine$integer.max - 2

 conBinFile <- file(description = binInvfile, open = "rb")
 Inv <- NULL
 while (length( a<-   readBin(conBinFile, 'double', MaxIntVal ) )> 0 ){
        Inv <- c(Inv, a)
}
 close(conBinFile)

 Inv <- matrix(data=Inv, nrow=nrow(Xmat) , byrow=FALSE)

    if(printInfo){
      message("Binary data has been read. Inverse matrix is now being formed. \n")
    }

  }



  return(Inv)

 }
