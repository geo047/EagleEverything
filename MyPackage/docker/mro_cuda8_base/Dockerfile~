# Author: Josh Bowden
# Company: CSIRO
# Date: 05/10/2017
# Description: NVidia Drivers 375.66 and CUDA SDK 8.0.44
# Image:  imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8

#### This docker file builds:
# docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base:latest

#### How did the image get there:
# sudo docker login imtsc-cont-reg.it.csiro.au
# sudo docker build -t mro_cuda8_base .
# sudo docker ps -a
# get the CONTAINER_ID for the image that was built (requires user to run the image first)
# export DPS=8aa27af3a63c
# sudo docker commit $DPS bow355/mro_cuda8_base
# sudo docker tag bow355/mro_cuda8_base  imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base
# sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base
# Retrive the image:
# sudo docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_base:latest


FROM imtsc-cont-reg.it.csiro.au/eagle/mro_shiny_base
# FROM mro_shiny_base_v3
MAINTAINER Josh Bowden "josh.bowden@csiro.au"

# R_LIBS_USER is defined in the base image
# ENV R_LIBS_USER /R/library
USER root

RUN  apt-get update \
    && apt-get install -y --no-install-recommends \  
        linux-headers-$(uname -r) \
        clinfo \
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/*
      
ENV DRIVERINSTALLDIR /usr/local/nvidia

# Install dependencies and Download and install shiny server
RUN  mkdir shinybuild \
    && LOCALIP=172.17.0.1 \
    && cd shinybuild \
    && NVID_VER=375.66 \ 
    && wget -O NVIDIA-Linux-x86_64-${NVID_VER}.run http://${LOCALIP}/nvidia/NVIDIA-Linux-x86_64-${NVID_VER}.run	\
    && chmod 755 NVIDIA-Linux-x86_64-${NVID_VER}.run \
    && ./NVIDIA-Linux-x86_64-${NVID_VER}.run --extract-only \
    && cd NVIDIA-Linux-x86_64-${NVID_VER} \
    && ln -s libEGL_nvidia.so.${NVID_VER}         libEGL_nvidia.so.0 \
    && ln -s libGLESv1_CM_nvidia.so.${NVID_VER}   libGLESv1_CM_nvidia.so.1 \
    && ln -s libGLESv2_nvidia.so.${NVID_VER}      libGLESv2_nvidia.so.2 \
    && ln -s libGLX_nvidia.so.${NVID_VER}         libGLX_indirect.so.0 \
    && ln -s libGLX_nvidia.so.${NVID_VER}         libGLX_nvidia.so.0 \
    && ln -s libnvidia-cfg.so.1                   libnvidia-cfg.so \
    && ln -s libnvidia-cfg.so.${NVID_VER}         libnvidia-cfg.so.1 \
    && ln -s libnvidia-encode.so.1                libnvidia-encode.so \
    && ln -s libnvidia-encode.so.${NVID_VER}      libnvidia-encode.so.1 \
    && ln -s libnvidia-fbc.so.1                   libnvidia-fbc.so \
    && ln -s libnvidia-fbc.so.${NVID_VER}         libnvidia-fbc.so.1 \
    && ln -s libnvidia-ifr.so.1                   libnvidia-ifr.so \
    && ln -s libnvidia-ifr.so.${NVID_VER}         libnvidia-ifr.so.1 \
    && ln -s libnvidia-ml.so.1                    libnvidia-ml.so \
    && ln -s libnvidia-ml.so.${NVID_VER}          libnvidia-ml.so.1 \
    && ln -s libnvidia-opencl.so.${NVID_VER}      libnvidia-opencl.so.1 \
    && ln -s vdpau/libvdpau_nvidia.so.${NVID_VER} libvdpau_nvidia.so \
    && ln -s libcuda.so.${NVID_VER}               libcuda.so \
    && ln -s libcuda.so.${NVID_VER}               libcuda.so.1 \
    && cd .. \
    && mkdir /usr/local/nvidia \
    && chmod -R 777 ${DRIVERINSTALLDIR} \
    && mv NVIDIA-Linux-x86_64-${NVID_VER}/* ${DRIVERINSTALLDIR} \
    # && ln -s /usr/local/nvidia/nvidia-smi  /usr/local/bin/nvidia-smi
    && mkdir -p /etc/OpenCL/vendors \
    && cp ${DRIVERINSTALLDIR}/nvidia.icd /etc/OpenCL/vendors	 \
    && find /etc/OpenCL \( -type d -exec chmod u+rwx,g+rwx,o+rwx {} \; -o -type f -exec chmod u+rw,g+rw,o+rw {} \; \) \
    && echo "${DRIVERINSTALLDIR}" >> /etc/ld.so.conf.d/cudadriver.conf \
    && ldconfig \
    && find ${DRIVERINSTALLDIR} \( -type d -exec chmod u+rwx,g+rwx,o+rwx {} \; -o -type f -exec chmod u+rw,g+rw,o+rw {} \; \) \
    && cd .. \
    && rm -rf shinybuild  


    # PATH includes /usr/local/cuda-8.0/bin
    # LD_LIBRARY_PATH includes /usr/local/cuda-8.0/lib64, 
    # /usr/local/cuda is created as a soft link to /usr/local/cuda-8.0/
    # chmod -R 777 /usr/local/cuda
RUN mkdir shinybuild \
    && cd shinybuild \
    && CUDA_VER=8.0.44 \
    LOCALIP=172.17.0.1 \
    && wget -O cuda_${CUDA_VER}_linux.run http://${LOCALIP}/nvidia/cuda_${CUDA_VER}_linux.run \
    && chmod 755 cuda_${CUDA_VER}_linux.run \
    && ./cuda_${CUDA_VER}_linux.run --silent --toolkit --no-opengl-libs --verbose --override \
    && echo "/usr/local/cuda-8.0/lib64" >> /etc/ld.so.conf \
    && ldconfig \
    && find /usr/local/cuda-8.0 \( -type d -exec chmod u+rwx,g+rwx,o+rwx {} \; -o -type f -exec chmod u+rw,g+rw,o+rw {} \; \) \
    && cd .. \
    && rm -rf shinybuild  
    

RUN  mkdir /opt/nvblas \
    && chmod 777 /opt/nvblas   

ADD nvblas.conf   /opt/nvblas/
ENV NVBLAS_CONFIG_FILE /opt/nvblas/nvblas.conf
RUN echo $WHEREISMAKECONF \
    && sed -i "s,__WHEREISMAKECONF__,$WHEREISMAKECONF,g" $NVBLAS_CONFIG_FILE \
    && chown docker:docker $NVBLAS_CONFIG_FILE \
    && chmod 666 $NVBLAS_CONFIG_FILE
    
    
# I will do this in the docker layer that will want to use the LD_PRELOAFD NV lib option
# ENV LD_PRELOAD /usr/local/cuda/lib64/libnvblas.so

ENV LD_LIBRARY_PATH="${DRIVERINSTALLDIR}:/usr/local/cuda-8.0/lib64:${LD_LIBRARY_PATH}"
ENV PATH="/usr/local/cuda-8.0/bin:${DRIVERINSTALLDIR}:${PATH}"

VOLUME  /data /flush1 /flush2

ENV    DATADIR /data 
ENV    FLUSHDIR /flush1 
ENV    FLUSH1DIR /flush1 
ENV    FLUSH2DIR /flush2 
ENV    MEMDIR /memdir 

USER docker

# add data files from local filesytem into the container
CMD ["R", "--vanilla"]


