.print_final  <- function(selected_loci, map,  extBIC, gamma )
{
  ## internal function: used by AM
  if (length(selected_loci) == 1 & any(is.na(selected_loci)))
  {
      message("No significant marker-trait associations have been found. \n\n")
     if (is.null(gamma)){
        message("Gamma value for model selection was set to its default value of 1. ")
     } else {
        message("Gamma value for model selection was set to  ", round(gamma,3) , " ")
     } 

  }  else {
     .print_results(selected_loci=selected_loci, map=map,  extBIC=extBIC)
     if (is.null(gamma)){
        message("Gamma value for model selection was set to its default value of 1. ")
     } else {
        message("Gamma value for model selection was set to  ", round(gamma,3) , " ")
     } 
  }   ## end if else
  message("")

}  ## end function print.finals


