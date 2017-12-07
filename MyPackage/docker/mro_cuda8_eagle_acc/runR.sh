#!/bin/bash
# Author: Josh Bowden
# Company: CSIRO
# Date: 07/12/2017
# Description: Script to launch R with or without LD_PRELOAD which is used to provide GPU BLAS Level 3 functions
# This script is part of the Eagle pakckage Docer version:
#  docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2:latest
# Usage: cat rscript.R | runR.sh <NUM_GPUS>

# if not being read from stdin then read a line
if  [ ! -t 0 ] ; then
    read PIPEDDATA
    echo "$PIPEDDATA"
fi 
# exit 0

if [ -z $1 ] ; then
  if [ ! -z "$PIPEDDATA" ] ; then
    echo "$PIPEDDATA" | R --vanilla
  else
    R --vanilla
  fi
    exit $?
fi

if [ "$1" -eq 0 ] ; then
   if [ ! -z "$PIPEDDATA" ] ; then
    echo "$PIPEDDATA" | R --vanilla
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
    if [ ! -z "$PIPEDDATA" ] ; then
     echo "$PIPEDDATA" | R --vanilla
    else
      R --vanilla
    fi
    exit $?
fi
