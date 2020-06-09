library(Eagle)

# Scenario 0
# Full data 
# 102 samples for pheno1
# Z matrix of 102 rows - first 3 records are repeats
# contains no missing trait data
# Turning age into factor and testing summary function 
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno1.test", missing="NA")
pheno$age <- as.factor(pheno$age)
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- AM(trait="trait1", fformula="age", Z=Z, map=map, geno=geno, pheno=pheno, lambda=0.5, 
fixit=TRUE, maxit=5)
SummaryAM(res)


Table 1: Summary Information 
   
--------------------------------------------------------
Number cpu:                               28        
Max memory (Gb):                          16        
Number of samples:                        102       
Number of snp:                            4998      
Trait name:                               trait1    
Fixed model:                              age                           
Number samples missing obs:               0         
Number significant snp-trait assocs:      4         
Gamma value for extBIC:                   0.5 
--------------------------------------------------------



 Number of cores being used for calculation is .. 28

 Table 2: Size and Significance of Effects in Final Model 
   
                       Effect Size     Df   Wald statstic       Pr(Chisq)
         (Intercept)         23.49      1           50.26       1.346E-12
               age19        -23.45      1           23.50       1.247E-06
               age20        -23.07      1           22.88       1.724E-06
               age21        -19.98      1           27.81       1.342E-07
               age22        -28.36      1           34.79       3.672E-09
               age23        -21.03      1           25.67       4.053E-07
               age24        -17.78      1           19.21       1.172E-05
               age25        -18.05      1           24.17       8.802E-07
               age26        -22.44      1           38.95       4.348E-10
               age27        -20.50      1           31.72       1.784E-08
               age28        -18.46      1           26.99       2.047E-07
               age29        -16.14      1           21.02       4.550E-06
               age30        -18.11      1           26.06       3.302E-07
               age31        -18.00      1           25.94       3.521E-07
               age32        -15.76      1           19.30       1.119E-05
               age33        -14.79      1           15.93       6.570E-05
               age34        -18.46      1           19.82       8.531E-06
               age35        -16.54      1           19.20       1.177E-05
               age37        -14.78      1           16.78       4.203E-05
               age38        -16.56      1           20.18       7.043E-06
               age40        -20.63      1           27.28       1.764E-07
               age41        -14.41      1           14.04       1.789E-04
               age42        -15.91      1           11.58       6.654E-04
               age55        -20.75      1           18.76       1.485E-05
           rs9974282          2.66      1           17.00       3.742E-05
         rs150917106          2.67      1           23.35       1.349E-06
           rs6122453          3.37      1           27.29       1.754E-07
           rs2079251          2.39      1           22.82       1.783E-06






# Scenario 1
# Full data 
# 102 samples for pheno1
# Z matrix of 102 rows - first 3 records are repeats
# contains no missing trait data
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno1.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)
SummaryAM(res)



 Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  0.995 
0.3157895  |  0.96 
0.3684211  |  0.88 
0.4210526  |  0.72 
0.4736842  |  0.575 
0.5263158  |  0.48 
0.5789474  |  0.4 
0.6315789  |  0.285 
0.6842105  |  0.225 
0.7368421  |  0.135 
0.7894737  |  0.085 
0.8421053  |  0.045 
0.8947368  |  0.03 
0.9473684  |  0.03 
1  |  0.02 
 ----------------------------- 
 For a false positive rate of  0.05  set the lambda parameter in the AM function to  0.8421053 

 New results after iteration 2 are 
            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         633.78   
     rs10404933          19          19035425           3899            636.79   
                           Final  Results  
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No significant marker-trait associations have been found. 


 Table 1: Summary Information 
   
--------------------------------------------------------
Number cpu:                               28        
Max memory (Gb):                          16        
Number of samples:                        102       
Number of snp:                            4998      
Trait name:                               trait1    
Fixed model:                              pc1 + pc2                     
Number samples missing obs:               0         
Number significant snp-trait assocs:      0         
Gamma value for extBIC:                   0.842105263157895
--------------------------------------------------------



 No significant marker-trait associations have been found by AM. 

 No p-values to report 






# Scenario 2
# Full data but with a missig trait value 
# 102 samples for pheno2 but there is one trait value missing and 2 covariate values missing.
# Data still contains repeat measures. 
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno2.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)


  Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  0.99 
0.3157895  |  0.955 
0.3684211  |  0.88 
0.4210526  |  0.73 
0.4736842  |  0.595 
0.5263158  |  0.48 
0.5789474  |  0.375 
0.6315789  |  0.295 
0.6842105  |  0.22 
0.7368421  |  0.145 
0.7894737  |  0.08 
0.8421053  |  0.04 
0.8947368  |  0.035 
0.9473684  |  0.03 
1  |  0.02 
 ----------------------------- 
 For a false positive rate of  0.05  set the lambda parameter in the AM function to  0.8421053 

 Significant marker-trait association found. 
 New results after iteration 2 are 

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         638.26   
     rs10404933          19          19035425           3899            640.98   

                           Final  Results  
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No significant marker-trait associations have been found. 



# Scenario 3
# 102 samples for pheno3 with records 7, 3, 2, 1 missing and record 9 missing covariate data
# Missing data causes pheno not to contain any repeat measures
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno3.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)

 Table: Empirical false positive rates, given lambda value for model selection. 

  Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  1 
0.3157895  |  0.96 
0.3684211  |  0.895 
0.4210526  |  0.745 
0.4736842  |  0.625 
0.5263158  |  0.51 
0.5789474  |  0.365 
0.6315789  |  0.255 
0.6842105  |  0.165 
0.7368421  |  0.12 
0.7894737  |  0.08 
0.8421053  |  0.04 
0.8947368  |  0.03 
0.9473684  |  0.025 
1  |  0.015 
 ----------------------------- 
 For a false positive rate of  0.05  set the lambda parameter in the AM function to  0.8421053 
New results after iteration 2 are 

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         651.45   
     rs10404933          19          19035425           3899            655.01   

                           Final  Results  
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No significant marker-trait associations have been found. 








# Scenario 4
# Results should be equivalent to above  results 
# 97 samples for pheno4 with no missing records
# no repeated measures. 
# 
library(Eagle)
geno <- ReadMarker(filename="./geno4.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno4.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2",  map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2",  map=map, geno=geno, pheno=pheno, lambda=res$setlambda)

#  Gamma    |  False Pos Rate  
# ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  1 
0.3157895  |  0.995 
0.3684211  |  0.965 
0.4210526  |  0.825 
0.4736842  |  0.695 
0.5263158  |  0.585 
0.5789474  |  0.475 
0.6315789  |  0.36 
0.6842105  |  0.255 
0.7368421  |  0.18 
0.7894737  |  0.105 
0.8421053  |  0.085 
0.8947368  |  0.045 
0.9473684  |  0.025 
1  |  0.025 
# ----------------------------- 
# For a false positive rate of  0.05  set the lambda parameter in the AM function to  0.8947368 
#
#           SNP        Chrm           Map Pos     Col Number            extBIC 
#          -----      ------         ---------     -----------         --------- 
#     Null Model                                                         607.50   
#     rs10404933          19          19035425           3899            612.61   







 


# Scenario 5
# 102 samples for pheno5 with  missing records
# trait data missing for 4 and 5. Covariate data also missing. 
# 
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)

 Table: Empirical false positive rates, given lambda value for model selection. 

  Gamma    |  False Pos Rate  
 ---------------------------- 





# Scenario 6
       library(Eagle)
       complete.name <- system.file('extdata', 'geno.ped', package='Eagle')
       geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
       complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
       pheno_obj <- ReadPheno(filename=complete.name)
       complete.name <- system.file('extdata', 'map.txt', package='Eagle')
       map_obj <- ReadMap(filename=complete.name)
       res <- AM(trait = 'y', fformula=c("cov1 + cov2"), map = map_obj, pheno = pheno_obj, geno = geno_obj)
       SummaryAM(AMobj=res)

 Number of cores being used for calculation is .. 28


 Table 1: Significance of Effects in Final Model 
   
                         Df   Wald statstic      Pr(Chisq)
         (Intercept)      1           3.20       7.367E-02
                cov1      1           0.92       3.365E-01
                cov2      1           0.83       3.610E-01
         rs145002694      1         120.57       0.000E+00
          rs77659166      1          49.84       1.665E-12




 Table 2: Proportion of phenotypic variance explained by the 
          model. Marker loci, which were found by AM(), are added
          a SNP at a time.
 
                   SNP          Proportion 
        +  rs145002694               0.402
         +  rs77659166               0.551


### Test cases with no fixed effects 

# Scenario 1
# Full data 
# 102 samples for pheno1
# Z matrix of 102 rows - first 3 records are repeats
# contains no missing trait data
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno1.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)

           SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         632.49   
      rs6122453          20          62826188           3230            639.35   




# Scenario 2
# Full data but with a missig trait value 
# 102 samples for pheno2 but there are 3 single records missing
# Data still contains repeat measures. 
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno2.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
    Null Model                                                         636.85   
      rs6122453          20          62826188           3230            643.94   







# Scenario 3
# Results  WILL NOT be equivalent to next lot of  results 
# 102 samples for pheno3 but records 1, 5, and 7 are missing
# Missing data causes pheno not to contain any repeat measures
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno3.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         650.61   
      rs6122453          20          62826188           3230            656.24   






# Scenario 4
# 97 samples for pheno4 with no missing records
# no repeated measures. 
# 
library(Eagle)
geno <- ReadMarker(filename="./geno4.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno4.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
res <- FPR4AM(trait="trait1",   map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1",   map=map, geno=geno, pheno=pheno, lambda=res$setlambda)


           SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         606.35   
      rs9974282          21          42333729           1189            612.99   
Gamma value for model selection was set to  0.842 






# Scenario 5
# 102 samples for pheno5 with  missing records
# all missing data [2,3,7,10,11] for single observation cases
#  repeated measures still present . 
# 
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, numreps=30)
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, lambda=res$setlambda)


 New results after iteration 2 are 
            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         637.49   
      rs1000728          21           3276368             90            645.16   

                           Final  Results  
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No significant marker-trait associations have been found. 
Gamma value for model selection was set to  0.842 





# Scenario 6
### Summary testings ... 
       library(Eagle)
       complete.name <- system.file('extdata', 'geno.ped', package='Eagle')
       geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8)
       complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
       pheno_obj <- ReadPheno(filename=complete.name)
       complete.name <- system.file('extdata', 'map.txt', package='Eagle')
       map_obj <- ReadMap(filename=complete.name)
       res <- AM(trait = 'y',  map = map_obj, pheno = pheno_obj, geno = geno_obj)
       SummaryAM(AMobj=res)


 Table 1: Summary Information 
   
--------------------------------------------------------
Number of samples:                        150       
Number of snp:                            100       
Trait name:                               y         
Fixed model:                              intercept only                
Number samples missing obs:               0         
Number significant snp-trait assocs:      2         

--------------------------------------------------------

 Number of cores being used for calculation is .. 28

 Table 2: Significance of Effects in Final Model 
   
                         Df   Wald statstic      Pr(Chisq)
         (Intercept)      1         486.17       0.000E+00
         rs145002694      1         118.29       0.000E+00
          rs77659166      1          48.22       3.813E-12




 false positive rate of  0.05  set the lambda parameter in the AM function to  0.7894737 
            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         638.61   
     rs10404933          19          19035425           3899            641.90   

No significant marker-trait associations have been found. 
Gamma value for model selection was set to  0.789 




# Scenario 5
# 102 samples for pheno5 with  missing records
# all missing data [2,3,7,10,11] for single observation cases
#  repeated measures still present . 
# lambda set to 0.3
# 
library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, lambda=0.36, maxit=10)
SummaryAM(AMobj=res)

                          Final  Results  

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         637.49   
      rs1000728          21           3276368             90            636.95   
      rs1350291          19          36305471           4508            631.29   
    rs113647042          20          40871298           2559            629.03   
    rs570109966          19           4140867           3382            626.70   
      rs6509094          19          44548677           4775            623.70   
     rs10404933          19          19035425           3899            621.03   
      rs6105077          20           1838137           1356            617.47   
     rs11699657          20          44582365           2675            608.32   
    rs199551210          21          37136022           1054            602.56   



Gamma value for model selection was set to  0.36 

 Table 1: Summary Information 
   
--------------------------------------------------------
Number cpu:                               28        
Max memory (Gb):                          16        
Number of samples:                        102       
Number of snp:                            4998      
Trait name:                               trait1    
Fixed model:                              intercept only                
Number samples missing obs:               2         
Number significant snp-trait assocs:      9         
Gamma value for extBIC:                   0.36
--------------------------------------------------------



 Number of cores being used for calculation is .. 28


 Table 2: Size and Significance of Effects in Final Model 
   
                       Effect Size     Df   Wald statstic       Pr(Chisq)
         (Intercept)          7.82      1           68.95       1.110E-16
           rs1000728          3.31      1           49.24       2.265E-12
           rs1350291         -2.66      1           25.05       5.598E-07
         rs113647042         -2.79      1           37.58       8.775E-10
         rs570109966          2.24      1           28.87       7.726E-08
           rs6509094         -2.84      1           28.93       7.522E-08
          rs10404933         -2.92      1           31.17       2.363E-08
           rs6105077         -2.27      1           28.74       8.295E-08
          rs11699657         -3.13      1           31.00       2.580E-08
         rs199551210          2.02      1           16.08       6.065E-05



