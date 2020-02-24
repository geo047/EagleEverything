#!/bin/bash


for f in *.cpp 
do 
echo "Processing $f file.." 
sed "s/asciifileM/tmpM/g" $f > tmp 
mv tmp $f
done
