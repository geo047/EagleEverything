#!/bin/bash

echo 

# cd /flush1/bow355/AMplus_new_code/Large
FILEPREFIX=eagle_HDF5_cs_eigenblas
FILEPREFIX=eagle_REPEATS_HDF5_ss
FILEPREFIX=eagle_REPEAT_cs_CRAN
FILEPREFIX=eagle_REPEAT_cs_ROWMAJOR
# list all files - possibly should not use the "-l" option here
# as then I would not need to remove the EOL in the section below
ls -l $FILEPREFIX*.res > ss.files
## Get a list of just the file names
cat ss.files | sed 's/.*eagle_/eagle_/' > ss2.files
## get a list of the SLURM_ID that is embeded in the file name
cat ss2.files | sed "s/${FILEPREFIX}_//" | sed "s/__.*//" > ss2_slurmid.files
# And a list of the individual repeat SLURM_ID
cat ss2.files | sed 's/.*__//' | sed 's/_gpu.*//' > ss2_slurmid_2.files

# This removes the line breaks so we can make an array
IFS=$'\r\n' GLOBIGNORE='*' command eval  'SLURMID_1=($(cat ss2_slurmid.files))'
IFS=$'\r\n' GLOBIGNORE='*' command eval  'SLURMID_2=($(cat ss2_slurmid_2.files))'
IFS=$'\r\n' GLOBIGNORE='*' command eval  'FILENAMES=($(cat ss2.files))'

# for NTHREADS in "${THREADLIST[@]}"; do
LOOPNUM=0
for FILES in "${FILENAMES[@]}"; do

    ADDID1="${SLURMID_1[${LOOPNUM}]}"
    ADDID2="${SLURMID_2[${LOOPNUM}]}"
    echo $ADDID1 $ADDID2
    
    cat $FILES | grep profile, | sed "s/profile,/profile,$ADDID1,$ADDID2,/g" |  sed "s/profile,$ADDID1,$ADDID2,itnum/profile,slurm_id1,slurm_id2,itnum/g" | sort -h -r > ${FILES%.*}.sid
    
    LOOPNUM=$(( $LOOPNUM + 1 ))
done

# then we require: 
REMOVELINES=`cat ${FILEPREFIX}_*.sid | grep unction,time_ms | wc -l` 
cat ${FILEPREFIX}_*.sid | sort -h -r | tail -n +${REMOVELINES} > ${FILEPREFIX}_slurm_id.txt
echo "Your results data is named: ${FILEPREFIX}_slurm_id.txt"  
echo "This script is not needed as Eagle now stores a run_id with the SLURM_JOB_ID"