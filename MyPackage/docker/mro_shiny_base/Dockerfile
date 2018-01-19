# Author: Josh Bowden
# Company: CSIRO
# Date: 05/10/2017

#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_shiny_base:latest



# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
FROM ubuntu:16.04

## User creation code (based on https://github.com/rocker-org/rocker/blob/master/r-base/Dockerfile):
## Set a default user. Available via runtime flag `--user docker`
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (e.g. for linked volumes to work properly).
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get clean && apt-get update && apt-get install -y locales \
    && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Install some useful tools and dependencies for MRO
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ca-certificates \
		curl \
		nano \
		# MRO dependencies dpkg does not install on its own:
		libcairo2 \
		libgfortran3 \
		libglib2.0-0 \
		libgomp1 \
		libjpeg8 \
		libpango-1.0-0 \
		libpangocairo-1.0-0 \
		libtcl8.6 \
		libtcl8.6 \
		libtiff5 \
		libtk8.6 \
		libx11-6 \
		libxt6 \
		# needed for installation of MKL:
		build-essential \
		make \
		gcc \
        apt-utils \
        apt-transport-https \
        sudo \
        pandoc \
        pandoc-citeproc \
        libcurl4-gnutls-dev \
        libcairo2-dev \
        libxt-dev \
        g++ \
        gfortran \
        bash-completion \
        libjpeg-dev \
        wget \
        git \
        libssl-dev \
        libssh2-1-dev \
        libssl1.0.0 \
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/*


## https://mran.revolutionanalytics.com/documents/rro/installation/#revorinst-lin
# Use major and minor vars to re-use them in non-interactive installation script
ENV MRO_VERSION_MAJOR 3
ENV MRO_VERSION_MINOR 4.2
ENV MRO_VERSION $MRO_VERSION_MAJOR.$MRO_VERSION_MINOR
# N.B. Change the .libPaths() on line 168 below to match the MRO version '/opt/microsoft/ropen/3.4.2/lib64/R/library'

WORKDIR /home/docker
ENV    WHEREISMAKECONF=/opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc

# Download & Install MRO
RUN    LOCALIP=172.17.0.1 \
    && wget -O  microsoft-r-open-$MRO_VERSION.tar.gz http://$LOCALIP/mro/microsoft-r-open-$MRO_VERSION.tar.gz \
    && tar -xzf microsoft-r-open-$MRO_VERSION.tar.gz \
    && cd microsoft-r-open/ \
    && ./install.sh -u -a \
    && cd .. \
    && rm -rf microsoft-r-open \
    && rm -rf microsoft-r-open-$MRO_VERSION.tar.gz
 
# RUN   ls -l /usr/lib64 \
#       && sleep 10 

# Print MKL license on every start 
# && RUN echo 'cat("\n", readLines("/home/docker/microsoft-r-open/MKL_EULA.txt"), "\n", sep="\n")' >> $WHEREISMAKECONF/Rprofile.site  
# WORKDIR /home/docker 


ARG VCS_URL
ARG VCS_REF
ARG BUILD_DATE

# Add image metadata
LABEL org.label-schema.license="https://mran.microsoft.com/faq/#licensing" \
    org.label-schema.vendor="Microsoft Corporation, based on Dockerfile of Daniel Nüst (https://hub.docker.com/r/nuest/mro/~/dockerfile/) and modified by Josh Bowden" \
	org.label-schema.name="Microsoft R Open and R Shiny packages" \
	org.label-schema.description="Docker images of Microsoft R Open (MRO) with the Intel® Math Kernel Libraries (MKL) and Shiny packages installed in /R/library." \ 
	org.label-schema.vcs-url=$VCS_URL \
	org.label-schema.vcs-ref=$VCS_REF \
	org.label-schema.build-date=$BUILD_DATE \
	org.label-schema.schema-version="rc1" \
	maintainer="Josh Bowden <josh.bowden@csiro.au>"


RUN   mkdir -p /data /shinydata /flush1 /flush2 /memdir /R /R/library /datastore /tmp/eagle_git /home/docker/.R \
   && chmod 777 /R/library 
   
# N.B. Change the .libPaths() on line 168 below to match the MRO version '/opt/microsoft/ropen/3.4.2/lib64/R/library'
# Remove this command in the Rprofile.site file as it casues an error in Docker image (not in Singularity though)
 RUN sed -i 's/utils::assignInNamespace("quit"/# utils::assignInNamespace("quit"/g' $WHEREISMAKECONF/Rprofile.site \
     && echo "suppressPackageStartupMessages(suppressWarnings(library(tcltk)))" >> $WHEREISMAKECONF/Rprofile.site \
     && echo ".libPaths( c( '/R/library', '/opt/microsoft/ropen/3.4.2/lib64/R/library') )" >> $WHEREISMAKECONF/Rprofile.site \
     && echo "local({ options(shiny.port = 3838, shiny.host = \"0.0.0.0\") })" >> $WHEREISMAKECONF/Rprofile.site 

     
USER docker 

ENV  MAKE='make -j2'
RUN  R  -e "install.packages(c('devtools', 'httpuv', 'shiny', 'shinythemes', 'shinyBS', 'shinyjs', 'rmarkdown') , lib='/R/library' , dependencies=TRUE)" 
RUN  R  -e "install.packages(c('matrixcalc', 'data.table', 'RcppEigen', 'Rcpp'),  lib='/R/library', dependencies=TRUE)" 


ENV   R_LIBS_USER /R/library

# /usr/bin/R -> /usr/lib64/microsoft-r/3.4/lib64/R/bin/R
# CMD  ["R", "-e", "library(Eagle); OpenGUI()"]
# R -e "library(Eagle); OpenGUI()"
# ENTRYPOINT  ["R", "-e", "library(Eagle); OpenGUI()"]
CMD  ["R", "--vanilla"]
