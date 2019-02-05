.form_results <- function(trait, selected_loci,   fformula, indxNA_pheno, indxNA_geno,
                           ncpu, availmemGb, quiet,  extBIC, gamma, geno, pheno, map, Zmat )
{
  ## internal function - used by AM for forming the results object
  if (length(selected_loci) > 1){
   sigres <- list(trait=trait,
                    fformula = fformula,
                    indxNA_pheno = indxNA_pheno,
                    indxNA_geno = indxNA_geno,
                    Mrk=map[[1]][selected_loci],
                    Chr=map[[2]][selected_loci],
                    Pos=map[[3]][selected_loci],
                    Indx=selected_loci,
                    ncpu=ncpu,
                    availmemGb=availmemGb,
                    quiet=quiet,
                    extBIC=extBIC,
                    gamma=gamma,
                    geno=geno,
                    pheno=pheno,
                    map=map,
                    Zmat=Zmat)
  } else {
   sigres <- list(trait=trait,
                    fformula = fformula,
                    indxNA_pheno = indxNA_pheno,
                    indxNA_geno = indxNA_geno,
                    Mrk=NA,
                    Chr=NA,
                    Pos=NA,
                    Indx=selected_loci,
                    ncpu=ncpu,
                    availmemGb=availmemGb,
                    quiet=quiet,
                    extBIC=extBIC,
                    gamma=gamma,
                    geno=geno,
                    pheno=pheno,
                    map=map,
                    Zmat=Zmat)
  }
return(sigres)
}


