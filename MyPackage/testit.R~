library(Eagle)

geno <- ReadMarker(filename="./genoDemo.dat", type="text", AA=0, AB=1, BB=2)

pheno <- ReadPheno(filename = "./phenoDemo.dat")

map <- ReadMap(filename = "./mapDemo.dat")


res <- FPR4AM(trait="trait1", fformula="pc1+pc2", map=map, geno=geno, pheno=pheno, falseposrate=0.05, numreps=250)


res <- AM(trait="trait1", fformula="pc1+pc2", map=map, geno=geno, pheno=pheno, gamma=res$setgamma)





#
#                           Final  Results
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#            SNP        Chrm           Map Pos     Col Number            extBIC
#          -----      ------         ---------     -----------         ---------
#     Null Model                                                         928.87
#      rs6132876          20          25857279           2207            884.29
#     rs12975028          19          35180369           4503            880.77
#
#
#
#Gamma value for model selection was set to  0.895
#



