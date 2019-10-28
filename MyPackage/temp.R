library(Eagle)
      complete.name <- system.file('extdata', 'map.txt', 
                                        package='Eagle')
       map_obj <- ReadMap(filename=complete.name) 
     
       # read marker data
       #~~~~~~~~~~~~~~~~~~~~
       # Reading in a PLINK ped file 
       # and setting the available memory on the machine for the reading of the data to 8  gigabytes
       complete.name <- system.file('extdata', 'geno.ped', 
                                          package='Eagle')
       geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
      
       # read phenotype data
       #~~~~~~~~~~~~~~~~~~~~~~~
     
       # Read in a plain text file with data on a single trait and two covariates
       # The first row of the text file contains the column names y, cov1, and cov2. 
       complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
       
       pheno_obj <- ReadPheno(filename=complete.name)
                
     
      #  Suppose we want to perform the AM analysis at a 5% false positive rate. 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
       ans <- FPR4AM(falseposrate = 0.05,
                     trait = 'y',
                     numreps=500,
                     fformula=c('cov1+cov2'),
                     map = map_obj,
                     pheno = pheno_obj,
                     geno = geno_obj) 
      
      res <- AM(trait =  'y',
                     fformula=c('cov1+cov2'),
                     map = map_obj,
                     pheno = pheno_obj,
                     geno = geno_obj,
                     gamma = ans$setgamma)
 

