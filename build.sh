#!/bin/bash

# build this project
mkdir build
cd build
# make clean
cmake -DMKL_ROOT_PATH=/home/george/intel/ -DMETIS_ROOT_PATH=/home/george/intel-5.1.0/build/Linux-x86_64/  -DCMAKE_BUILD_TYPE=Release ..
make
