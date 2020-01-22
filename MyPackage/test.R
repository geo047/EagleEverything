library(Lion)


geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2, availmemGb=0.005 )
pheno <- ReadPheno(filename = "./pheno1.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- FPR4AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno)
