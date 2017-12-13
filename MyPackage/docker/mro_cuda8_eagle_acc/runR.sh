#!/bin/bash
# Author: Josh Bowden
# Company: CSIRO
# Date: 07/12/2017
# Description: Script to launch R with or without LD_PRELOAD which is used to provide GPU BLAS Level 3 functions
# This script is part of the Eagle pakckage Docer version:
#  docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2:latest
# Usage: cat rscript.R | runR.sh [NUM_GPUS]
# N.B. the rscript.R file must have a carriage return at last line
# cat am.R | EAGLE_PROFILE_STR=1  OMP_NUM_THREADS=14  singularity  exec --nv mro_cuda8_eagle_acc2.img /usr/bin/runR.sh

# if not being read from stdin then read in all lines and save to variable RFILE
if  [ ! -t 0 ] ; then
    RFILE=""
    while read -r -t 10 PIPEDDATA;
    do
      RFILE+="${PIPEDDATA}\n"
    done
fi

# echo -e "$RFILE"
# exit 

export EAGLE_PROFILE_STR=1 
export OMP_NUM_THREADS=$OMP_NUM_THREADS


echo "We will be using the following version of R:"
which R
# exit 0

if [ -z $1 ] ; then
  if [ ! -z "$RFILE" ] ; then
    echo -e "$RFILE" | R --vanilla
  else
    R --vanilla
  fi
    exit $?
fi

if [ "$1" -eq 0 ] ; then
   if [ ! -z "$RFILE" ] ; then
    echo -e "$RFILE" | R --vanilla
  else
    R --vanilla
  fi
    exit $?
fi

GPUSTOUSE=`seq -s "," 0 $(( $1 - 1 ))` 
# GPUSTOUSE=`seq -s " " 0 $(( 4 - 1 ))` 
if [ "$1" -ge 1 ] ; then
    # sed -i "s/NVBLAS_GPU_LIST 0/NVBLAS_GPU_LIST $GPUSTOUSE/g" $NVBLAS_CONFIG_FILE
    export CUDA_VISIBLE_DEVICES=$GPUSTOUSE
    export LD_PRELOAD=/usr/local/cuda/lib64/libnvblas.so
    #export NVBLAS_CONFIG_FILE=/flush1/bow355/AMplus_new_code/Mid_docker_tests/nvblas_01234567.conf
    #export LD_PRELOAD=$CUDA_HOME/lib64/libnvblas.so
    if [ ! -z "$RFILE" ] ; then
     echo -e "$RFILE" | R --vanilla
    else
      R --vanilla
    fi
    exit $?
fi

