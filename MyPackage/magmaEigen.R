magmaEigen <- function(Xmat , ngpu=1, wantvectors=TRUE, printInfo=FALSE){


## R frontend to Magma multi-GPU code for eigendecomposition for square symmetric matrix.
# Due to the 32 bit interface problem between R and Magma, the Magma code sits outside of R. It is 
# built separate to R. It is the executable that is being accessed here via a system() command. 
# The executable is in  the Magma directory.  

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
  maxnumrows <- sqrt( .Machine$integer.max/sizeof(double)) - 8 # for a bit of a buffer
  if ( nrow(Xmat) > maxnumrows ){  
    # write file in chunks by blocking the Xmat based on column blocks  
    # Step 1: determine how many columns will fit in MaxIntVal
    numcols <- trunc( MaxIntVal/( nrow(Xmat)*sizeof(double) ) )
    if (numcols < 0){
       message(" Error: number of rows of data matrix is larger than the max integer value of the machine. ")
       return( 0)
    }

    numblocks <- trunc( ncol(Xmat) / numcols)   
    extracols <- 0
    if (  (  ncol(Xmat) / numcols ) %% 1 > 0)
     {
           # this are extra values to be read in 
           extracols <- (  (  ncol(Xmat) /numcols  ) %% 1) * numcols
     }


     conBinFile <- file(description = binXmatfile, open = "wb")
      if (printInfo){
           message(" The Xmat matrix is large and must be written out  as ", numblocks+1, " blocks \n")
      }
      for(ii in 1:numblocks){
        indx <- seq( (ii-1)*numcols+1, ii*numcols)
        cat("Block number ", ii, " min indx ", min(indx), " max indx " , max(indx), "\n")
        if(printInfo){
           message(" Writing out block number ", ii, " of ", numblocks+1, " blocks of binary data. \n")
        }
     writeBin(object=as.vector(Xmat[, indx]), con=conBinFile)
      }
      if(extracols > 0){
        if(printInfo){
           message(" Writing out block number ", ii+1, " of ", numblocks+1, " blocks of binary data. \n")
        }
        indx <- seq(numblocks*numcols + 1, ncol(Xmat))
        cat(" Extra - min indx ", min(indx), " max indx " , max(indx), "\n")
        writeBin(object=as.vector(Xmat[, indx]), con=conBinFile)
     }


  } else  {
     writeBin(object=as.vector(Xmat), con=binXmatfile)
  }


  complete.name <- system.file('Magma', 'magma_eigen.exe', package='EagleGPU')
  if (printInfo)
    print(complete.name)
  system(paste(complete.name, binXmatfile, nrow(Xmat), ngpu, as.numeric(printInfo), binvalfile, binvecfile, as.numeric(wantvectors), "&> output.out" ) ) 

   # read in eigenvalues
   values <- readBin(binvalfile , double(), nrow(Xmat) )
   
   if (wantvectors){ 

     MaxIntVal <- .Machine$integer.max - 2

     conBinFile <- file(description = binvecfile, open = "rb")
     vectors <- NULL
     while (length( a<-   readBin(conBinFile, 'double', MaxIntVal ) )> 0 ){
          vectors <- c(vectors, a)
     }
     close(conBinFile)

     vectors <- matrix(data=vectors, nrow=nrow(Xmat) , byrow=FALSE)


   }

 
  # need to change the order of the results to be consistent with R order
  indx <- order(values, decreasing=TRUE)


  return(list(vectors=vectors[,indx], values=values[indx]))

 }



