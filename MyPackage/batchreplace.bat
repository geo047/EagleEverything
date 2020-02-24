#!/bin/bash


for f in *.R 
do 
echo "Processing $f file.." 
sed "s/asciifileM/tmpM/g" $f > tmp 
mv tmp $f
done
