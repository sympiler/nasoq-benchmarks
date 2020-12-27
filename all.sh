#!/bin/bash

DATASET=SMP_Repository/
BUILDIR=build/
METIS_PATH=/home/zjming1/metis-5.1.0/

if [ "$#" -ge 3 ]; then
DATASET=$1
BUILDIR=$2
METIS_PATH=$3
fi
eps=-3
echo "Running solvers in $BUILDIR for QP problems in $DATASET ..."

# make a directory for building project
mkdir -p build
cd build

# build NASOQ
if [ ! -d nasoq ]; then
    cmake -DCMAKE_PREFIX_PATH=$MKLROOT/lib/intel64/ ..
    cmake -DCMAKE_PREFIX_PATH=$MKLROOT/include/ ..
    cmake -DCMAKE_PREFIX_PATH=${METIS_PATH}build/Linux-x86_64/libmetis/ ..
    cmake -DMETIS_ROOT_PATH=$METIS_PATH ..
fi

for eps in {-3,-6,-9}; do
    echo "Running NASOQ-Fixed ..."
    bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps> logs/nasoq-fixed-e${eps}.csv

    echo "Running NASOQ-Tuned ..."
    bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v tuned"> logs/nasoq-tuned-e${eps}.csv

    echo "Running customized NASOQ ..."
    bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v predet -r 0"> logs/nasoq-custom-e${eps}.csv
done