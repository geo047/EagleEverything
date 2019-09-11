magmaQR <- function(Xmat , ngpu=1, printInfo=FALSE){

## an Magma function to calculate the Q matrix from QR factorisatino 
## of the square matrix X

if (nrows(X) != ncols(X)){
   message("\n magmaQR function needs a square matrix. \n")
   return(0)
}

# A temporary file name is assigned for the cpp function to use
 binQfile <- paste(tempdir() , "/", "Q.bin", sep="")


 succes <- magma_qr(X=Xmat ,  numgpus=ngpu, printInfo=printInfo, fname=binQfile )
 if(success !=0{
   messge("\n magmaQR function has failed with error ", success, "\n")
   return(0)
  } else {

  Q <- readBin(binQfile , double(), nrows(X) * ncols(X) )
  Q <- matrix(data=Q, nrow=nrows(X) , byrow=FALSE)
  }


  return(Qmat)

}
