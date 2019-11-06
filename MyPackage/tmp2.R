library(Eagle)
geno <- ReadMarker(filename="./geno.test", type="text", AA=0, AB=1, BB=2 )
pheno <- ReadPheno(filename = "./pheno2.test", missing="NA")
map <- ReadMap(filename = "./mapDemo.dat")
Z <- ReadZmat("./Z1.test")
res <- AM(trait="trait1", fformula="pc1+pc2", Z=Z, map=map, geno=geno, pheno=pheno, gamma=1, quiet=FALSE)



