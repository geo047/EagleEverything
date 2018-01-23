# Author: Josh Bowden
# Company: CSIRO
# Date: 05/10/2017

#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_master:latest

#### How did the image get there:
# sudo docker login imtsc-cont-reg.it.csiro.au
# sudo docker build -t mro_eagle_master .
# sudo docker ps -a
# get the CONTAINER_ID for the image that was built (requires user to run the image first)
# export DPS=8aa27af3a63c
# sudo docker commit $DPS bow355/mro_eagle_master
# sudo docker tag bow355/mro_eagle_master  imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_master
# sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_master
# Retrive the image:
# sudo docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_master:latest


FROM imtsc-cont-reg.it.csiro.au/eagle/mro_shiny_base
MAINTAINER Josh Bowden "josh.bowden@csiro.au"

# R_LIBS_USER is defined in the base image
# ENV R_LIBS_USER /R/library
USER root

ENV EAGLEDIR /R/library/Eagle/shiny_app

# RUN  R  -e "install.packages(c('shinyFiles'),  lib='/R/library', dependencies=TRUE)"  # this is done in mro_shiny_base image

# Install dependencies and Download and install shiny server
RUN  mkdir shinybuild \
    && cd shinybuild \
    &&  git clone https://github.com/geo047/EagleEverything . \
    &&  cd MyPackage   \
    &&  NPROCS=`getconf _NPROCESSORS_ONLN`  \
    &&  MAKE="make -j$NPROCS" R CMD INSTALL Eagle -l /R/library  \
    # R -e "install.packages(c('Eagle'), repos='http://cran.rstudio.com/' , lib='/R/library')" 
    &&  cd ../..  \
    &&  rm -rf shinybuild  \
    &&  mkdir -p $EAGLEDIR\shinydata \
    &&  echo "local({ options(shiny.port = 3838, shiny.host = \"0.0.0.0\") })" >> $WHEREISMAKECONF/Rprofile.site  \
    &&  echo "library(shinyFiles, lib.loc=\"/R/library\")" >> $WHEREISMAKECONF/Rprofile.site 
    # &&  echo "export R_LIBS_USER=/R/library" >> /home/docker/.bashrc 


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
ADD shinydata/*.dat   $EAGLEDIR/shinydata/

CMD ["R", "-e shiny::runApp('/R/library/Eagle/shiny_app')"]

