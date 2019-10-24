.print_final  <- function(selected_loci, map,  extBIC, gamma )
{
  ## internal function: used by AM
  if (length(selected_loci) == 1 & any(is.na(selected_loci)))
  {
      message("No significant marker-trait associations have been found. \n\n")
  }  else {
     .print_results(selected_loci=selected_loci, map=map,  extBIC=extBIC)
     if (is.null(gamma)){
        message("Gamma value for model selection was set to its default value of 1. \n")
     } else {
        message("Gamma value for model selection was set to  ", round(gamma,3) , " \n")
     } 
     message("\n\n")
  }   ## end if else


}  ## end function print.finals


