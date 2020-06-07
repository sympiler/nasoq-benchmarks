#!/bin/bash

TOOLBIN=$1
DATAPATH=$2

echo "Running $TOOLBIN for dataset in $DATAPATH ..."

QPP=$(find $DATAPATH -name "*.yml"  -type f)


for f in $QPP; do
	$TOOLBIN -i $f ;
done 