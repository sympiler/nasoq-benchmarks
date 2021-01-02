#!/bin/bash

DATASET=test_smp/
MKL_PATH=$MKLROOT/
METIS_PATH=/home/zjming1/metis-5.1.0/

if [ "$#" -ge 3 ]; then
DATASET=$1
MKL_PATH=$2
METIS_PATH=$3
fi
eps=-3
echo "Running solvers in $BUILDIR for QP problems in $DATASET ..."

# make a directory for building project
if [ -d build ]; then
    rm -rf build
fi
mkdir -p build
cd build

# build a makefile
if [ ! -f Makefile ]; then
    cmake -DCMAKE_PREFIX_PATH=$MKL_PATH/lib/intel64/ ..
    cmake -DCMAKE_PREFIX_PATH=$MKL_PATH/include/ ..
    cmake -DCMAKE_PREFIX_PATH=${METIS_PATH}build/Linux-x86_64/libmetis/ ..
    cmake -DMETIS_ROOT_PATH=$METIS_PATH ..
fi

# build the project
if [ -d nasoq ]  && [ -d drivers ]; then
    make
fi

# change to the root directory
cd ..

for eps in {-3,-6}; do
    echo "Running NASOQ-Fixed ..."
    bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $DATASET $eps> logs/nasoq-fixed-e${eps}.csv

    echo "Running NASOQ-Tuned ..."
    bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $DATASET $eps "-v tuned"> logs/nasoq-tuned-e${eps}.csv

    echo "Running customized NASOQ ..."
    bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $DATASET $eps "-v predet -r 0"> logs/nasoq-custom-e${eps}.csv

    echo "Running OSQP ..."
    bash scripts/NASOQ_bench.sh build/drivers/osqp-bench $DATASET $eps> logs/osqp-e${eps}.csv

    echo "Running OSQP-polished ..."
    bash scripts/NASOQ_bench.sh build/drivers/osqp-bench $DATASET $eps "-v polished"> logs/osqp-polished-e${eps}.csv
done

if [ python -c 'import sys; print(sys.version_info[0])' -eq 2 ]; then
    python3 all.py logs
else
    python all.py logs
fi