library(Lion)

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
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)
SummaryAM(res)



 Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  1 
0.3157895  |  0.98 
0.3684211  |  0.93 
0.4210526  |  0.81 
0.4736842  |  0.705 
0.5263158  |  0.59 
0.5789474  |  0.475 
0.6315789  |  0.34 
0.6842105  |  0.245 
0.7368421  |  0.16 
0.7894737  |  0.105 
0.8421053  |  0.05 
0.8947368  |  0.03 
0.9473684  |  0.03 
1  |  0.02 
 ----------------------------- 
 For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8421053 

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
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno2.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)


  Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  1 
0.3157895  |  0.98 
0.3684211  |  0.93 
0.4210526  |  0.82 
0.4736842  |  0.69 
0.5263158  |  0.595 
0.5789474  |  0.455 
0.6315789  |  0.355 
0.6842105  |  0.24 
0.7368421  |  0.17 
0.7894737  |  0.105 
0.8421053  |  0.045 
0.8947368  |  0.035 
0.9473684  |  0.03 
1  |  0.02 
 ----------------------------- 
 For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8421053 

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
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno3.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

 Table: Empirical false positive rates, given gamma value for model selection. 

  Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  1 
0.3157895  |  0.975 
0.3684211  |  0.915 
0.4210526  |  0.83 
0.4736842  |  0.715 
0.5263158  |  0.575 
0.5789474  |  0.42 
0.6315789  |  0.315 
0.6842105  |  0.18 
0.7368421  |  0.135 
0.7894737  |  0.095 
0.8421053  |  0.05 
0.8947368  |  0.03 
0.9473684  |  0.025 
1  |  0.015 
 ----------------------------- 
 For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8421053 
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
library(Lion)
geno <- ReadMarker(filename="./geno4.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno4.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2",  map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2",  map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

#  Gamma    |  False Pos Rate  
# ---------------------------- 
#0  |  1 
#0.05263158  |  1 
#0.1052632  |  1 
#0.1578947  |  1 
#0.2105263  |  1 
#0.2631579  |  1 
#0.3157895  |  0.995 
#0.3684211  |  0.98 
#0.4210526  |  0.935 
#0.4736842  |  0.83 
#0.5263158  |  0.695 
#0.5789474  |  0.55 
#0.6315789  |  0.42 
#0.6842105  |  0.31 
#0.7368421  |  0.225 
#0.7894737  |  0.14 
#0.8421053  |  0.105 
#0.8947368  |  0.055 
#0.9473684  |  0.035 
#1  |  0.025 
# ----------------------------- 
# For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8947368 
#
#           SNP        Chrm           Map Pos     Col Number            extBIC 
#          -----      ------         ---------     -----------         --------- 
#     Null Model                                                         607.50   
#     rs10404933          19          19035425           3899            612.61   







 


# Scenario 5
# 102 samples for pheno5 with  missing records
# trait data missing for 4 and 5. Covariate data also missing. 
# 
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

 Table: Empirical false positive rates, given gamma value for model selection. 

  Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  1 
0.2105263  |  1 
0.2631579  |  1 
0.3157895  |  0.985 
0.3684211  |  0.97 
0.4210526  |  0.9 
0.4736842  |  0.73 
0.5263158  |  0.59 
0.5789474  |  0.44 
0.6315789  |  0.32 
0.6842105  |  0.25 
0.7368421  |  0.13 
0.7894737  |  0.075 
0.8421053  |  0.05 
0.8947368  |  0.04 
0.9473684  |  0.02 
1  |  0.02 
 ----------------------------- 
 For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8421053 
 New results after iteration 2 are 
            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         638.61   
     rs10404933          19          19035425           3899            642.79   

                           Final  Results  
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No significant marker-trait associations have been found. 






# Scenario 6
### Summary testings ... 
       library(Lion)
       complete.name <- system.file('extdata', 'geno.ped', package='Lion')
       geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
       complete.name <- system.file('extdata', 'pheno.txt', package='Lion')
       pheno_obj <- ReadPheno(filename=complete.name)
       complete.name <- system.file('extdata', 'map.txt', package='Lion')
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
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

           SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         632.49   
      rs6122453          20          62826188           3230            639.35   




# Scenario 2
# Full data but with a missig trait value 
# 102 samples for pheno2 but there are 3 single records missing
# Data still contains repeat measures. 
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno2.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
    Null Model                                                         636.85   
      rs6122453          20          62826188           3230            643.94   







# Scenario 3
# Results  WILL NOT be equivalent to next lot of  results 
# 102 samples for pheno3 but records 1, 5, and 7 are missing
# Missing data causes pheno not to contain any repeat measures
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno3.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         650.61   
      rs6122453          20          62826188           3230            656.24   






# Scenario 4
# 97 samples for pheno4 with no missing records
# no repeated measures. 
# 
library(Lion)
geno <- ReadMarker(filename="./geno4.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno4.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
res <- FPR4AM(trait="trait1",   map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1",   map=map, geno=geno, pheno=pheno, gamma=res$setgamma)


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
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, numreps=30)
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)


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
       library(Lion)
       complete.name <- system.file('extdata', 'geno.ped', package='Lion')
       geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8)
       complete.name <- system.file('extdata', 'pheno.txt', package='Lion')
       pheno_obj <- ReadPheno(filename=complete.name)
       complete.name <- system.file('extdata', 'map.txt', package='Lion')
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








# Scenario 5
# 102 samples for pheno5 with  missing records
# all missing data [2,3,7,10,11] for single observation cases
#  repeated measures still present . 
# gamma set to 0.3
# 
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, gamma=0.36)
SummaryAM(AMobj=res)



 Table 1: Summary Information 
   
--------------------------------------------------------
Number of samples:                        102       
Number of snp:                            4998      
Trait name:                               trait1    
Fixed model:                              intercept only                
Number samples missing obs:               2         
Number significant snp-trait assocs:      20        
Gamma value for extBIC:                   0.36      
--------------------------------------------------------


 Table 2: Significance of Effects in Final Model 
   
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
     rs62119935          19          20518272           3954            597.01   
    rs113645990          19          48787304           4937            588.85   
     rs13040618          20          11758325           1660            577.42   
     rs71352280          19          46597270           4846            569.25   
      rs6054947          20           9481072           1587            566.91   
       rs478723          19           5234328           3427            562.17   
      rs7409488          19           9254491           3572            558.43   
      rs9979212          21          20985093            615            555.07   
    rs145676387          21          20530009            605            551.29   
    rs150917106          19          11082925           3621            547.22   


