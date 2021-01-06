#!/bin/bash

# load all modules
module load cmake/3.12.3
module load gcc/7.3.0
module load intel/2019.3
module load python/3.7.4
module load scipy-stack/2019b

DATASET=SMP_Repository/
MKL_PATH=$MKLROOT
METIS_PATH=/home/zjming1/metis-5.1.0/

echo "$#"
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
    cmake -DCMAKE_PREFIX_PATH=${MKL_PATH}/lib/intel64/ -DCMAKE_BUILD_TYPE=Release  ..
    cmake -DCMAKE_PREFIX_PATH=${MKL_PATH}/include/ -DCMAKE_BUILD_TYPE=Release  ..
    cmake -DCMAKE_PREFIX_PATH=${METIS_PATH}build/Linux-x86_64/libmetis/ -DCMAKE_BUILD_TYPE=Release ..
    cmake -DMETIS_ROOT_PATH=${METIS_PATH} -DCMAKE_BUILD_TYPE=Release ..
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

var="$(python -c 'import sys; print(sys.version_info[0])')"
if [[ $var == 2 ]]; then
    python3 all.py logs
else
    python all.py logs
fi