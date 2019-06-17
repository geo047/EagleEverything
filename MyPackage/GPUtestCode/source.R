# this needs to run in R


library(Rcpp)  


PKG_LIBS <- sprintf('-Wl,-rpath, /apps/magma/2.5.1a1-cuda90/lib -L /apps/magma/2.5.1a1-cuda90/lib -lmagma  /apps/magma/2.5.1a1-cuda90/lib/libmagma.a -Wl,-rpath,/apps/cuda/8.0.61/lib64 %s', Rcpp:::RcppLdFlags()) 
PKG_CPPFLAGS <- sprintf('-DADD_ -DHAVE_CUBLAS -I/apps/magma/2.5.1a1-cuda90//include -I/apps/cuda/8.0.61/include %s', Rcpp:::RcppCxxFlags())  


#PKG_CPPFLAGS <- sprintf(' -D_MAGMA_WITH_GPUS=1 -DADD_ -DMAGMA_WITH_MKL -openmp -fpic  -I/apps/magma/2.5.1a1-cuda90/include -I/apps/magma/2.5.1a1-cuda90/control -DHAVE_CUBLAS -I/apps/cuda/8.0.61/include  -I/home/geo047/RLibs/Rcpp/include/')
#PKG_LIBS <- sprintf('  /apps/R/3.4.0/lib64/R/bin/Rscript   -fPIC -openmp -Wl,-rpath,/apps/magma/2.5.1a1-cuda90/lib/lp64:/apps/intel/mkl/2017.2.174/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64 -L/apps/magma/2.5.1a1-cuda90/lib/lp64 -L/apps/cuda/8.0.61/lib64 -lpthread -lmagma -lcudart -lcublas -lcuda -lstdc++ -lm  ')

print("PKG_LIBS")
print(PKG_LIBS)

print("")
print("PKG_CPPFLAGS")
print(PKG_CPPFLAGS)


#PKG_LIBS <- sprintf('-Wl,-rpath,/usr/local/magma/lib -L/usr/local/magma/lib -lmagma /usr/local/magma/lib/libmagma.a -Wl,-rpath,/usr/local/cuda-5.5/lib64 %s', Rcpp:::RcppLdFlags()) 
#PKG_CPPFLAGS <- sprintf('-DADD_ -DHAVE_CUBLAS -I$(CUDA_HOME)/include  -I$(MAGMA_HOME)/include  %s', Rcpp:::RcppCxxFlags())  
Sys.setenv(PKG_LIBS = PKG_LIBS , PKG_CPPFLAGS = PKG_CPPFLAGS) 
R <- file.path(R.home(component = 'bin'), 'R') 
file <- '~/Test/gpuQR_magma.cpp'
cmd <- sprintf('%s CMD SHLIB %s', R, paste(file, collapse = ' '))
system(cmd)


## After compilation of the shared library ... 

dyn.load('./gpuQR_magma.so')
set.seed(100)
n_row <- 3; n_col <- 3
A <- matrix(rnorm(n_row * n_col), n_row, n_col)
qr(A)$qr
.Call('gpuQR_magma', A)



