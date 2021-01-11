#!/bin/bash

# load all modules
module load cmake
module load gcc
module load intel
module load python
module load scipy-stack

DATASET=SMP_Repository/
MKL_PATH=$MKLROOT
METIS_PATH=/home/zjming1/metis-5.1.0/

if [ "$#" -ge 3 ]; then
DATASET=$1
MKL_PATH=$2
METIS_PATH=$3
fi
eps=-3
echo "Running solvers in $BUILDIR for QP problems in $DATASET ..."

source ${MKL_PATH}/bin/mklvars.sh intel64

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
    cmake -DCMAKE_PREFIX_PATH=${METIS_PATH}/build/Linux-x86_64/libmetis/ -DCMAKE_BUILD_TYPE=Release ..
    cmake -DMETIS_ROOT_PATH=${METIS_PATH}/ -DCMAKE_BUILD_TYPE=Release ..
fi

# build the project
if [ -d nasoq ]  && [ -d drivers ]; then
    make
fi

# change to the root directory
cd ..

# make a directory for performance data
if [ -d logs/perf_data ]; then
    rm -rf logs/perf_data
fi
mkdir -p logs/perf_data
cd logs/perf_data

for eps in {-3,-6}; do
    for d in ${DATASET}/*/; do
        echo "Running NASOQ-Fixed ..."
        bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $d $eps> logs/perf_data/nasoq-fixed-e${eps}-${d}.csv

        echo "Running NASOQ-Tuned ..."
        bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $d $eps "-v tuned"> logs/perf_data/nasoq-tuned-e${eps}-${d}.csv

        echo "Running customized NASOQ ..."
        bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $d $eps "-v predet -r 0"> logs/perf_data/nasoq-custom-e${eps}-${d}.csv

        echo "Running OSQP ..."
        bash scripts/NASOQ_bench.sh build/drivers/osqp-bench $d $eps> logs/perf_data/osqp-e${eps}-${d}.csv

        echo "Running OSQP-polished ..."
        bash scripts/NASOQ_bench.sh build/drivers/osqp-bench $d $eps "-v polished"> logs/perf_data/osqp-polished-e${eps}-${d}.csv
    done

    for f in ${DATASET}/*.csv; do
        echo "Running NASOQ-Fixed ..."
        bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $f $eps> logs/perf_data/nasoq-fixed-e${eps}-non-class.csv

        echo "Running NASOQ-Tuned ..."
        bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $f $eps "-v tuned"> logs/perf_data/nasoq-tuned-e${eps}-non-class.csv

        echo "Running customized NASOQ ..."
        bash scripts/NASOQ_bench.sh build/nasoq/NASOQ-BIN $f $eps "-v predet -r 0"> logs/perf_data/nasoq-custom-e${eps}-non-class.csv

        echo "Running OSQP ..."
        bash scripts/NASOQ_bench.sh build/drivers/osqp-bench $f $eps> logs/perf_data/osqp-e${eps}-non-class.csv

        echo "Running OSQP-polished ..."
        bash scripts/NASOQ_bench.sh build/drivers/osqp-bench $f $eps "-v polished"> logs/perf_data/osqp-polished-e${eps}-non-class.csv
    done
done

var="$(python -c 'import sys; print(sys.version_info[0])')"
if [[ $var == 2 ]]; then
    python3 all.py logs/perf_data
else
    python all.py logs/perf_data
fi