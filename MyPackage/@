# LD_PRELOAD=$CUDA_ROOT/lib64/libnvblas.so R

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

 Table: Empirical false positive rates, given gamma value for model selection. 

#  Gamma    |  False Pos Rate  
# ---------------------------- 
#0  |  1 
#0.05263158  |  1 
#0.1052632  |  1 
#0.1578947  |  1 
#0.2105263  |  1 
#0.2631579  |  1 
#0.3157895  |  0.98 
#0.3684211  |  0.93 
#0.4210526  |  0.81 
#0.4736842  |  0.705 
#0.5263158  |  0.59 
#0.5789474  |  0.475 
#0.6315789  |  0.34 
#0.6842105  |  0.245 
#0.7368421  |  0.16 
#0.7894737  |  0.105 
#0.8421053  |  0.05 
#0.8947368  |  0.03 
#0.9473684  |  0.03 
#1  |  0.02 
# ----------------------------- 
# For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8421053 
# New results after iteration 2 are 
#
#            SNP        Chrm           Map Pos     Col Number            extBIC 
#          -----      ------         ---------     -----------         --------- 
#     Null Model                                                         633.78   
#     rs10404933          19          19035425           3899            636.79   







# Scenario 2
# Full data but with a missig trait value 
# 102 samples for pheno2 but there are 3 single records missing
# Data still contains repeat measures. 
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno2.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

#  Gamma    |  False Pos Rate  
# ---------------------------- 
#0  |  1 
#0.05263158  |  1 
#0.1052632  |  1 
#0.1578947  |  1 
#0.2105263  |  0.995 
#0.2631579  |  0.99 
#0.3157895  |  0.95 
#0.3684211  |  0.905 
#0.4210526  |  0.79 
#0.4736842  |  0.665 
#0.5263158  |  0.465 
#0.5789474  |  0.335 
#0.6315789  |  0.275 
#0.6842105  |  0.195 
#0.7368421  |  0.135 
#0.7894737  |  0.095 
#0.8421053  |  0.065 
#0.8947368  |  0.05 
#0.9473684  |  0.04 
#1  |  0.02 
# ----------------------------- 
# For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8947368 
#
# New results after iteration 2 are 
#
#            SNP        Chrm           Map Pos     Col Number            extBIC 
#          -----      ------         ---------     -----------         --------- 
#     Null Model                                                         618.00   
#     rs10404933          19          19035425           3899            622.35   
#






# Scenario 3
# Results should be equivalent to next lot of  results 
# 102 samples for pheno3 but records 1, 5, and 7 are missing
# Missing data causes pheno not to contain any repeat measures
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno3.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

#  Gamma    |  False Pos Rate  
# ---------------------------- 
#0  |  1 
#0.05263158  |  1 
#0.1052632  |  1 
#0.1578947  |  0.99 
#0.2105263  |  0.985 
#0.2631579  |  0.975 
#0.3157895  |  0.925 
#0.3684211  |  0.855 
#0.4210526  |  0.78 
#0.4736842  |  0.645 
#0.5263158  |  0.525 
#0.5789474  |  0.425 
#0.6315789  |  0.33 
#0.6842105  |  0.255 
#0.7368421  |  0.19 
#0.7894737  |  0.12 
#0.8421053  |  0.09 
#0.8947368  |  0.05 
#0.9473684  |  0.03 
#1  |  0.02 
# ----------------------------- 
# For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8947368 
#New results after iteration 2 are 
#
#            SNP        Chrm           Map Pos     Col Number            extBIC 
#          -----      ------         ---------     -----------         --------- 
#     Null Model                                                         607.50   
#     rs10404933          19          19035425           3899            612.61   







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
# all missing data [2,3,7,10,11] for single observation cases
#  repeated measures still present . 
# 
library(Lion)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno5.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)

 Gamma    |  False Pos Rate  
 ---------------------------- 
0  |  1 
0.05263158  |  1 
0.1052632  |  1 
0.1578947  |  0.995 
0.2105263  |  0.98 
0.2631579  |  0.945 
0.3157895  |  0.9 
0.3684211  |  0.84 
0.4210526  |  0.71 
0.4736842  |  0.6 
0.5263158  |  0.485 
0.5789474  |  0.355 
0.6315789  |  0.28 
0.6842105  |  0.19 
0.7368421  |  0.125 
0.7894737  |  0.075 
0.8421053  |  0.05 
0.8947368  |  0.02 
0.9473684  |  0.015 
1  |  0.005 
 ----------------------------- 
 For a false positive rate of  0.05  set the gamma parameter in the AM function to  0.8421053 

 New results after iteration 2 are 

            SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         602.95   
     rs10404933          19          19035425           3899            608.03   

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
     Null Model                                                         627.09   
      rs6122453          20          62826188           3230            634.29   






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
     Null Model                                                         611.52   
      rs6122453          20          62826188           3230            618.50   






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
      rs9974282          21          42333729           1189            613.88   






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
res <- FPR4AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno)
res <- AM(trait="trait1",  Z=Z, map=map, geno=geno, pheno=pheno, gamma=res$setgamma)


          SNP        Chrm           Map Pos     Col Number            extBIC 
          -----      ------         ---------     -----------         --------- 
     Null Model                                                         618.04   
      rs1000728          21           3276368             90            625.91   





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


 Table 1: Significance of Effects in Final Model 
   
                         Df   Wald statstic      Pr(Chisq)
         (Intercept)      1         486.17       0.000E+00
         rs145002694      1         118.29       0.000E+00
          rs77659166      1          48.22       3.813E-12




 Table 2: Proportion of phenotypic variance explained by the 
          model. Marker loci, which were found by AM(), are added
          a SNP at a time.
 
                   SNP          Proportion 
        +  rs145002694               0.403
         +  rs77659166               0.548

