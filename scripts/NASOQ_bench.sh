#!/bin/bash

TOOLBIN=$1
DATAPATH=$2
ACC=$3
header=1
#echo "Running $TOOLBIN for dataset in $DATAPATH ..."

QPP=$(find $DATAPATH -name "*.yml"  -type f)


for f in $QPP; do
 if [ $header -eq 1 ]; then
  $TOOLBIN -i $f -d 1 -e $ACC
  header=0
 else
  $TOOLBIN -i $f -e $ACC
 fi
 echo ""

done 
