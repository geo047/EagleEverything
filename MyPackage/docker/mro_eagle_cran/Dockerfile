# Author: Josh Bowden
# Company: CSIRO
# Date: 30/11/2017


#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_cran:latest

#### How did the image get there:
# sudo docker login imtsc-cont-reg.it.csiro.au
# sudo docker build -t mro_eagle_cran .
# sudo docker ps -a
# get the CONTAINER_ID for the image that was built (requires user to run the image first)
# export DPS=8aa27af3a63c
# sudo docker commit $DPS bow355/mro_eagle_cran
# sudo docker tag bow355/mro_eagle_cran  imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_cran
# sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_cran
# Retrive the image:
# sudo docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_cran:latest

FROM imtsc-cont-reg.it.csiro.au/eagle/mro_shiny_base
MAINTAINER Josh Bowden "josh.bowden@csiro.au"

# R_LIBS_USER is defined in the base image
# ENV R_LIBS_USER /R/library
USER root

ENV EAGLEDIR /R/library/Eagle/shiny_app
 

# Install dependencies and Download and install shiny server
RUN R -e "install.packages(c('Eagle'), lib='/R/library')" \
    &&  mkdir -p $EAGLEDIR\shinydata \
    &&  echo "library(Eagle)" >> $WHEREISMAKECONF/Rprofile.site \ 
    &&  echo "local({ options(shiny.port = 3838, shiny.host = \"0.0.0.0\") })" >> $WHEREISMAKECONF/Rprofile.site 
    # &&  echo "export R_LIBS_USER=/R/library" >> /home/docker/.bashrc 


EXPOSE 3838

USER docker

WORKDIR /R/library/Eagle/shiny_app

VOLUME  /data /flush1 /flush1 /flush2

ENV    DATADIR /data 
ENV    FLUSHDIR /flush1 
ENV    FLUSH1DIR /flush1 
ENV    FLUSH2DIR /flush2 
ENV    MEMDIR /memdir 

USER docker

# add data files from local filesytem into the container
ADD shinydata/*.dat   $EAGLEDIR/shinydata/

CMD ["R", "-e Eagle::OpenGUI()"]  

