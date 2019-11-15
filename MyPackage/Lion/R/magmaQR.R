magmaQR <- function(Xmat ,  printInfo=FALSE){


## The R frontend to the magma code (magma_qr) 
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
 binQfile <- paste(tempdir() , "/", "Q.bin", sep="")


 success <- magma_qr(X=Xmat ,  numgpus=computer$ngpu, printInfo=printInfo, fname=binQfile  )
 if(success < 0)
  {
   message("\n magmaQR function has failed. A zero value has been returned.  \n")
   return(0)
  } else {

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
      Q <- NULL
      conBinFile <- file(description = binQfile, open = "rb")
      if (printInfo){
           message(" The input matrix is large and must be read in as ", numblocks+1, " blocks \n")
      }
      for(ii in 1:numblocks){
        if(printInfo){
           message(" Reading in block number ", ii, " of ", numblocks+1, " blocks of binary data. \n") 
        }
        Q <- c(Q, readBin(conBinFile , double(), MaxIntVal))
      }
      if(extravals > 0){
        if(printInfo){
           message(" Reading in block number ", ii+1, " of ", numblocks+1, " blocks of binary data. \n") 
        }
        Q <-  c(Q, readBin(conBinFile , double(), extravals) )
     }

    if(printInfo){
      message("Binary data has been read. Q matrix is now being formed. \n")
    }


    Q <- matrix(data=Q, nrow=nrow(Xmat) , byrow=FALSE)

   } else {
      Q <- readBin(binQfile , double(), nrow(Xmat) * ncol(Xmat) )
      Q <- matrix(data=Q, nrow=nrow(Xmat) , byrow=FALSE)
   }


  }


  return(Q)

 }
