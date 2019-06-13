# Docker files and use of images for the Eagle package

Author: Josh Bowden
Company: CSIRO
Date: 06/12/2017

### CSIRO IMT Docker registry  

Available: [https://imtsc-cont-reg.it.csiro.au/harbor/eagle](https://imtsc-cont-reg.it.csiro.au/harbor/projects/21/repository)

The docker images have the follwing hierachy:

* mro_shiny_base       
   * mro_eagle_cran
   * mro_eagle_master
   * mro_eagle_acc2
   * mro_cuda8_base            
        * mro_cuda8_eagle_master
        * mro_cuda8_eagle_acc2
 
 N.B. 
 mro_shiny_base image requires MRO to be available for download via wgetat:
      http://$LOCALIP/mro/microsoft-r-open-$MRO_VERSION.tar.gz
 mro_cuda8_base image requires the GPU driver installation file and CUDA installation files need to be available at:
      http://${LOCALIP}/nvidia/NVIDIA-Linux-x86_64-${NVID_VER}.run
      http://${LOCALIP}/nvidia/cuda_${CUDA_VER}_linux.run
 
 
### To get the images from the registry
```
docker login imtsc-cont-reg.it.csiro.au
docker pull imtsc-cont-reg.it.csiro.au/eagle/mro_cuda8_eagle_acc2:latest
```

### To update the registry with new images
```
 sudo docker login imtsc-cont-reg.it.csiro.au
 # build an image
 sudo docker build --no-cache -t mro_shiny_base .
 sudo docker tag bow355/mro_eagle_cran  imtsc-cont-reg.it.csiro.au/eagle/mro_shiny_base
 sudo docker push imtsc-cont-reg.it.csiro.au/eagle/mro_shiny_base
```

### Running docker on Linux 
Pass through the DISPAY variable to the container so the tcltk file.chooser works
```
 sudo docker run --rm  -p 3838:3838 -it  mro_eagle_master
```
 
###  Use singularity to launch the container
```
module load singularity

SINGULARITY_DOCKER_USERNAME=<ident>
SINGULARITY_DOCKER_PASSWORD=
singularity build mro_eagle_master.simg docker://imtsc-cont-reg.it.csiro.au/eagle/mro_eagle_master

# export the mount points for the CSIRO clusters:
export SINGULARITY_BINDPATH="$DATADIR:/data,$FLUSHDIR:/flush1,$FLUSH1DIR:/flush1,$FLUSH2DIR:/flush2,$MEMDIR:/memdir"
 
# run the image then launch a web browser pointing to: http://localhost:3838
./mro_eagle_master.simg
 N.B. As singularity binds a users home directory, the shiny app will not find the 
 /home/docker/shinydata directory, so manual entry of the path is needed to run the embedded demo data model

# Or use the image as an R interpreter
# pipe an R script into the containers R process:
cat am+GPU.R | OMP_NUM_THREADS=14 singularity exec --nv mro_eagle_master.simg /usr/bin/runR.sh 0 --vanilla

# N.B. the --nv option should bind the GPU drivers into the container, however there are current issues
# with this when the drivers are in a non-standard place
# i.e. on the CSIRO clusters they are in /cm/local/apps/cuda-driver/libs/current/lib64
```
