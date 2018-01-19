# Author: Josh Bowden
# Company: CSIRO
# Date: 05/12/2017

#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_master:latest

#### How did the image get there:
# sudo docker login imtsc-cont-reg.it.csiro.au
# sudo docker build -t mro_cuda8_master .
# sudo docker ps -a
# get the CONTAINER_ID for the image that was built (requires user to run the image first)
# export DPS=8aa27af3a63c
# sudo docker commit $DPS bow355/mro_cuda8_master
# sudo docker tag bow355/mro_cuda8_master  imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_master
# sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_master
# Retrive the image:
# sudo docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_master:latest


FROM imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base
MAINTAINER Josh Bowden "josh.bowden@csiro.au"

# R_LIBS_USER is defined in the base image
# ENV R_LIBS_USER /R/library
USER root

ENV EAGLEDIR /R/library/Eagle/shiny_app


# Install dependencies and Download and install shiny server
RUN  mkdir shinybuild \
    && cd shinybuild \
    &&  git clone https://github.com/geo047/EagleEverything . \
    &&  cd MyPackage   \
    &&  cat Eagle/inst/shiny_app/app.R \
    &&  NPROCS=`getconf _NPROCESSORS_ONLN`  \
    &&  MAKE="make -j$NPROCS" R CMD INSTALL Eagle -l /R/library  \
    # R -e "install.packages(c('Eagle'), lib='/R/library')" 
    &&  cd ../..  \
    &&  rm -rf shinybuild  \
    &&  mkdir -p $EAGLEDIR\shinydata 
    # &&  echo "local({ options(shiny.port = 3838, shiny.host = \"0.0.0.0\") })" >> $WHEREISMAKECONF/Rprofile.site \
    # &&  echo "export R_LIBS_USER=/R/library" >> /home/docker/.bashrc 


EXPOSE 3838

RUN chmod 666 $NVBLAS_CONFIG_FILE
# ENV LD_PRELOAD /usr/local/cuda/lib64/libnvblas.so

USER docker

WORKDIR /R/library/Eagle/shiny_app


# add data files from local filesytem into the container
ADD shinydata/*.dat   $EAGLEDIR/shinydata/


CMD ["R", "-e Eagle::OpenGUI()"]  
