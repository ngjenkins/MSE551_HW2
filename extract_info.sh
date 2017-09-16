#!/bin/bash

declare N=25;
declare filename=1-HW2_1_;

cd $PWD

for i in `seq 0 $N`;
do
#################################################extract_info###########################################################

E=`grep  -A1 "Energy initial, next-to-last, final =" ${filename}_out/${filename}$i.out` 
echo $i $E  >>SUMMARY

#################################################end_extract_info#######################################################
done



