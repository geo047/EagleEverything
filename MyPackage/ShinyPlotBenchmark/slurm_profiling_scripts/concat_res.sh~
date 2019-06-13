#!/bin/bash

# FILEPREFIX=eagle_REPEAT_cs_ROWMAJOR
FILEPREFIX=$1
if [[ -z $FILEPREFIX ]] ; then
echo "Please enter the filename prefix on the command line"
exit -1
fi
REMOVELINES=`cat ${FILEPREFIX}_*.res | grep unction,time_ms | wc -l` 
cat ${FILEPREFIX}*.res | grep profile, | sort -r | tail -n +${REMOVELINES} > ${FILEPREFIX}_run_id.txt
echo "Your results data is named: ${FILEPREFIX}_run_id.txt" 