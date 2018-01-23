# Author: Josh Bowden
# Company: CSIRO
# Date: 05/12/2017
# Built images are available CSIRO Docker image repo: 

#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2:latest

#### How did the image get there:
# sudo docker login imtsc-cont-reg.it.csiro.au
# sudo docker build -t mro_cuda8_eagle_acc2_hdf .
# sudo docker tag mro_cuda8_eagle_acc2_hdf  imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2_hdf
# sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2_hdf
# Get image on cluster - https://confluence.csiro.au/display/eagleaccel/Run+Eagle+using+Singularity+on+a+Cluster
# 
# disable the cache as it is stored in your home directory which has limited space
# export SINGULARITY_DISABLE_CACHE=yes
# give your repository credentials (only needed if the image is not public)
# export SINGULARITY_DOCKER_USERNAME=<ident>
# Yes this is not the best thing to be doing
# export SINGULARITY_DOCKER_PASSWORD=<psswd>
# singularity pull docker://imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2_hdf:latest

# sudo docker ps -a
# get the CONTAINER_ID for the image that was built (requires user to run the image first)
# export DPS=8aa27af3a63c
# sudo docker commit $DPS bow355/mro_cuda8_eagle_acc2

FROM imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base
MAINTAINER Josh Bowden "josh.bowden@csiro.au"

# R_LIBS_USER is defined in the base image
USER root

ENV EAGLEDIR /R/library/Eagle/shiny_app

# RUN apt-get update \
#    && apt-get install -y --no-install-recommends \
#        libhdf5-serial-dev \
#     && apt-get clean \
#     && rm -rf /var/lib/apt/lists/*
#
     
#     # R -e "install.packages(c('Eagle'), repos='http://cran.rstudio.com/' , lib='/R/library')" 
# Install dependencies and Download and install shiny server
RUN  mkdir shinybuild \
    && cd shinybuild \
    && git clone https://github.com/geo047/EagleEverything . \
    && git checkout acc2  \
    && cd MyPackage/Eagle/src   \
    && sed -i 's/#include <RcppEigen.h>/#define EIGEN_USE_BLAS\n#include <RcppEigen.h>/g' *.cpp \
    && cd ../..  \
    && NPROCS=`getconf _NPROCESSORS_ONLN`  \
    && MAKE="make -j$NPROCS" R CMD INSTALL Eagle -l /R/library  \
    && cd ../..  \
    && rm -rf shinybuild  \
    && pwd \
    && mkdir -p $EAGLEDIR/shinydata 


EXPOSE 3838
ADD runR.sh /usr/bin

ADD nvblas_01234567.conf   /opt/nvblas/
ENV NVBLAS_CONFIG_FILE  /opt/nvblas/nvblas_01234567.conf 
RUN echo $WHEREISMAKECONF \
    && sed -i "s,__WHEREISMAKECONF__,$WHEREISMAKECONF,g" $NVBLAS_CONFIG_FILE \
    && chown docker:docker $NVBLAS_CONFIG_FILE \
    && chmod 666 $NVBLAS_CONFIG_FILE
    
RUN chmod 666 $NVBLAS_CONFIG_FILE \
    && chmod 777 /usr/bin/runR.sh
    
    
# ENV LD_PRELOAD /usr/local/cuda/lib64/libnvblas.so

USER docker

WORKDIR $EAGLEDIR/shinydata

VOLUME  /data /flush1 /flush1 /flush2

# add data files from local filesytem into the container
ADD shinydata/*.dat $EAGLEDIR/shinydata/


CMD ["R", "-e Eagle::OpenGUI()"] 
