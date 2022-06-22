#!/bin/bash                                                                                                     

CWD=$(pwd)

headernames=''
rm -f rms.out
numfiles=$(ls -l ????_*_*_*_13.pdb | wc -l)
for i in ????_*_*_*_13.pdb;  do 

    for j in ????_*_*_*_13.pdb;  do 
	b=$(basename "$i"); 
	/TCRmodeller/programs/ProFitV3.1/src/profit -f profit.in -h ${i} ${j} > profit.out
	var=$(grep RMS: profit.out | awk '{print $NF}')
	values=$values','$var 
    done
    headernames=""$headernames"",""$b""
done
echo "B = matrix( c($values), nrow=$numfiles,ncol=$numfiles)" >> rms.out
#echo -e $values > rms.out
echo "dimnames(A) = list( c($headernames), c($headernames)" >> rms.out
