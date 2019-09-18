magmaEigen <- function(Xmat , ngpu=1, wantvectors=TRUE, printInfo=FALSE){


## The R frontend to the magma code (magma_eigen) 
## Due to the 32 bit interface problem between R and Magma, the Q matrix solution from 
## magma_qr is passed back into R via a temporary binary file using the readBin() function. 
## An issue I encountered is that R only handles 32 bit integers, even for 64 bit builds. This 
## means large matrices are a problem for readBin.  To get around the problem, large Q matrices 
## are read  back into R as  blocks.  Fortunately, the reading of binary data is very fast. 
## If an error is encountered in magma_qr, most likely due to not enough memory, then this function 
## returns with a 0 value. 

 if (nrow(Xmat) != ncol(Xmat)){
   message("\n magmaQR function needs a square matrix. \n")
   return(0)
 }

 # A temporary file name is assigned for the cpp function to use
 binvalfile <-  paste(tempdir() , "/", "values.bin", sep="")

 binvecfile <- NULL
 if(wantvectors){
    binvecfile <-  paste(tempdir() , "/", "vectors.bin", sep="")
}


 # Arg - have to write the Xmat to disc to get this to work
  binXmatfile <-  paste(tempdir() , "/", "Xmat.bin", sep="")
  # will need to change this if this works for 32bit version
  writeBin(object=as.vector(Xmat), con=binXmatfile)

  complete.name <- system.file('Magma', 'magma_eigen.exe', package='Eagle')
  system(paste(complete.name, binXmatfile, nrow(Xmat), ngpu, as.numeric(printInfo), binvalfile, binvecfile, as.numeric(wantvectors)))



print("Success")  

# success <- magma_eigen(X=binXmatfile  , numrows=nrow(Xmat),  numgpus=ngpu, printInfo=printInfo, fnameval=binvalfile, fnamevec=binvecfile, message=message, wantvectors=wantvectors )
 if(success < 0)
  {
   message("\n magmaEigen function has failed. A zero value has been returned.  \n")
   return(0)
  } else {
   # read in eigenvalues
   values <- readBin(binvalfile , double(), nrow(Xmat) )
   values <- matrix(data=values, nrow=nrow(Xmat) , byrow=FALSE)
   vectors <- NULL 
   if (wantvectors){ 

     lognumvals <- log(nrow(Xmat)) + log(ncol(Xmat)) 
     MaxIntVal <- .Machine$integer.max - 8  # just for a little bit of safety

     if (  lognumvals  > log(MaxIntVal)  ) {
         numblocks <- trunc( ( nrow(Xmat) / MaxIntVal) * ncol(Xmat) )  # getting around 32 bit integer issues
         extravals <- 0
         if (  ( ( nrow(Xmat) / MaxIntVal) * ncol(Xmat) ) %% 1 > 0)
         {
           # this are extra values to be read in 
           extravals <- (  ( ( nrow(Xmat) / MaxIntVal) * ncol(Xmat) ) %% 1) * MaxIntVal

         } 
      vectors <- NULL
      conBinFile <- file(description = binvecfile, open = "rb")
      if (printInfo){
           message(" The input matrix is large and must be read in as ", numblocks+1, " blocks \n")
      }
      for(ii in 1:numblocks){
        if(printInfo){
           message(" Reading in block number ", ii, " of ", numblocks+1, " blocks of binary data. \n") 
        }
        vectors <- c(vectors, readBin(conBinFile , double(), MaxIntVal))
      }
      if(extravals > 0){
        if(printInfo){
           message(" Reading in block number ", ii+1, " of ", numblocks+1, " blocks of binary data. \n") 
        }
        vectors <-  c(vectors, readBin(conBinFile , double(), extravals) )
     }

    if(printInfo){
      message("Binary data has been read. vector matrix is now being formed. \n")
    }


    vectors <- matrix(data=vectors, nrow=nrow(Xmat) , byrow=FALSE)

   } else {
      vectors <- readBin(binvecfile , double(), nrow(Xmat) * ncol(Xmat) )
      vectors <- matrix(data=vectors, nrow=nrow(Xmat) , byrow=FALSE)
   }

  }
  }  ## end if else
 
  # need to change the order of the results to be consistent with R order
  indx <- order(values, decreasing=TRUE)

  return(list(vectors=vectors[,indx], values=values[indx]))

 }
