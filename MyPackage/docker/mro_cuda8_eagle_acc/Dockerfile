# Author: Josh Bowden
# Company: CSIRO
# Date: 05/12/2017
# Built images are available CSIRO Docker image repo: 

#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_acc2:latest

#### How did the image get there:
# sudo docker login imtsc-cont-reg.it.csiro.au
# sudo docker build -t mro_cuda8_acc2 .
# sudo docker ps -a
# get the CONTAINER_ID for the image that was built (requires user to run the image first)
# export DPS=8aa27af3a63c
# sudo docker commit $DPS bow355/mro_cuda8_acc2
# sudo docker tag bow355/mro_cuda8_acc2  imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_acc2
# sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_acc2


FROM imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base
MAINTAINER Josh Bowden "josh.bowden@csiro.au"

# R_LIBS_USER is defined in the base image
USER root

ENV EAGLEDIR /R/library/Eagle/shiny_app

#     # R -e "install.packages(c('Eagle'), repos='http://cran.rstudio.com/' , lib='/R/library')" 
# Install dependencies and Download and install shiny server
RUN  mkdir shinybuild \
    && cd shinybuild \
    && git clone https://github.com/geo047/EagleEverything . \
    && git checkout acc2  \
    && cd MyPackage   \
    && NPROCS=`getconf _NPROCESSORS_ONLN`  \
    && MAKE="make -j$NPROCS" R CMD INSTALL Eagle -l /R/library  \
    && cd ../..  \
    && rm -rf shinybuild  \
    && pwd \
    && mkdir -p $EAGLEDIR/shinydata 


EXPOSE 3838

RUN chmod 666 $NVBLAS_CONFIG_FILE
# ENV LD_PRELOAD /usr/local/cuda/lib64/libnvblas.so

USER docker

WORKDIR $EAGLEDIR/shinydata

VOLUME  /data /flush1 /flush1 /flush2

# add data files from local filesytem into the container
ADD shinydata/*.dat $EAGLEDIR/shinydata/


CMD ["R", "-e Eagle::OpenGUI()"] 
