
# dynamic build - in line with R

icpc  -std=c++11 -DADD_ -I$MAGMA_HOME/include -I$CUDA_HOME/include -m64 -DMKL_ILP64   -I${MKLROOT}/include   -c magma_eigennonsym.cpp
icpc  -L$MAGMA_HOME/lib  -lmagma -L$CUDA_HOME/lib64 -lcublas -lcudart -lcusparse -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl    -o magma_eigennonsym.exe  magma_eigennonsym.o

