magmaQR <- function(Xmat , ngpu=1, printInfo=FALSE){

## an Magma function to calculate the Q matrix from QR factorisatino 
## of the square matrix X

if (nrow(Xmat) != ncol(Xmat)){
   message("\n magmaQR function needs a square matrix. \n")
   return(0)
}

# A temporary file name is assigned for the cpp function to use
 binQfile <- paste(tempdir() , "/", "Q.bin", sep="")


 success <- magma_qr(X=Xmat ,  numgpus=ngpu, printInfo=printInfo, fname=binQfile, message=message )
 if(success !=0)
  {
   messge("\n magmaQR function has failed with error ", success, "\n")
   return(0)
  } else {

  lognumvals <- log(nrow(Xmat)) + log(ncol(Xmat)) 
  MaxIntVal <- .Machine$integer.max - 8  # just for a little bit of safety
   MaxIntVal <- 20000*20000 - 15000    # testing 
  cat("lognumvals = ", lognumvals, "\n")
  cat("MatIntVal = ", MaxIntVal, "\n")

#Q <-  readBin(binQfile , double(),  as.double(nrow(Xmat) * ncol(Xmat)) )
#Q <- matrix(data=Q, nrow=nrow(Xmat) , byrow=FALSE)
#print("your kidding ... ")

#  return(Q)

  if (  lognumvals  > log(MaxIntVal)  ) {
      numblocks <- trunc( ( nrow(Xmat) / MaxIntVal) * ncol(Xmat) )  # getting around 32 bit integer issues
      extravals <- 0
      if (  ( ( nrow(Xmat) / MaxIntVal) * ncol(Xmat) ) %% 1 > 0)
      {
        # this are extra values to be read in 
        extravals <- (  ( ( nrow(Xmat) / MaxIntVal) * ncol(Xmat) ) %% 1) * MaxIntVal

      } 
      cat("Number of block reads needed ", numblocks, "\n")
      Q <- NULL
      conBinFile <- file(description = binQfile, open = "rb")
      for(ii in 1:numblocks){
        print(ii)
        #indx <- seq((ii-1)*MaxIntVal+1 , ii*MaxIntVal )
        #cat("min and max indx values", min(indx), " " , max(indx), "\n")
        #Q[indx] <- readBin(binQfile , double(), MaxIntVal)
        Q <- c(Q, readBin(conBinFile , double(), MaxIntVal))
      }
      if(extravals > 0){
        cat("extravals = ",  extravals, "\n")
        #indx <- seq( (ii-1)*MaxIntVal+1, as.double(nrow(Xmat)*ncol(Xmat)) )
        Q <-  c(Q, readBin(conBinFile , double(), extravals) )
     }
    print(" about to form Q matrix")
    Q <- matrix(data=Q, nrow=nrow(Xmat) , byrow=FALSE)

   } else {
      Q <- readBin(binQfile , double(), nrow(Xmat) * ncol(Xmat) )
      Q <- matrix(data=Q, nrow=nrow(Xmat) , byrow=FALSE)
   }


  }


  return(Q)

}
