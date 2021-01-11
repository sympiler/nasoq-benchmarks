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

dir=$PWD
cd $DATASET

for eps in {-3,-6}; do
    for d in *; do
        echo $d
        if [ -d "$d" ]; then
            echo "Running NASOQ-Fixed ..."
            bash ${dir}/scripts/NASOQ_bench.sh ${dir}/build/nasoq/NASOQ-BIN $d $eps> ${dir}/logs/perf_data/nasoq-fixed-e${eps}-${d}.csv

            echo "Running NASOQ-Tuned ..."
            bash ${dir}/scripts/NASOQ_bench.sh ${dir}/build/nasoq/NASOQ-BIN $d $eps "-v tuned"> ${dir}/logs/perf_data/nasoq-tuned-e${eps}-${d}.csv

            echo "Running customized NASOQ ..."
            bash ${dir}/scripts/NASOQ_bench.sh ${dir}/build/nasoq/NASOQ-BIN $d $eps "-v predet -r 0"> ${dir}/logs/perf_data/nasoq-custom-e${eps}-${d}.csv

            echo "Running OSQP ..."
            bash ${dir}/scripts/NASOQ_bench.sh ${dir}/build/drivers/osqp-bench $d $eps> ${dir}/logs/perf_data/osqp-e${eps}-${d}.csv

            echo "Running OSQP-polished ..."
            bash ${dir}/scripts/NASOQ_bench.sh ${dir}/build/drivers/osqp-bench $d $eps "-v polished"> ${dir}/logs/perf_data/osqp-polished-e${eps}-${d}.csv
        else
            if [ "${d: -4}" == ".yml" ]; then
                echo "Running NASOQ-Fixed ..."
                if [ ! -f ${dir}/logs/perf_data/nasoq-fixed-e${eps}-non-class.csv ]; then
                    ${dir}/build/nasoq/NASOQ-BIN -i $d -e $eps -d 1 > ${dir}/logs/perf_data/nasoq-fixed-e${eps}-non-class.csv
                else
                    ${dir}/build/nasoq/NASOQ-BIN -i $d -e $eps >> ${dir}/logs/perf_data/nasoq-fixed-e${eps}-non-class.csv
                fi

                echo "Running NASOQ-Tuned ..."
                if [ ! -f ${dir}/logs/perf_data/nasoq-tuned-e${eps}-non-class.csv ]; then
                    ${dir}/build/nasoq/NASOQ-BIN -i $d -e $eps -v tuned -d 1 > ${dir}/logs/perf_data/nasoq-tuned-e${eps}-non-class.csv
                else
                    ${dir}/build/nasoq/NASOQ-BIN -i $d -e $eps -v tuned >> ${dir}/logs/perf_data/nasoq-tuned-e${eps}-non-class.csv
                fi

                echo "Running customized NASOQ ..."
                if [ ! -f ${dir}/logs/perf_data/nasoq-custom-e${eps}-non-class.csv ]; then
                    ${dir}/build/nasoq/NASOQ-BIN -i $d -e $eps -v predet -r 0 -d 1 > ${dir}/logs/perf_data/nasoq-custom-e${eps}-non-class.csv
                else
                    ${dir}/build/nasoq/NASOQ-BIN -i $d -e $eps -v predet -r 0 >> ${dir}/logs/perf_data/nasoq-custom-e${eps}-non-class.csv
                fi

                echo "Running OSQP ..."
                if [ ! -f ${dir}/logs/perf_data/osqp-e${eps}-non-class.csv ]; then
                    ${dir}/build/drivers/osqp-bench -i $d -e $eps -d 1 > ${dir}/logs/perf_data/osqp-e${eps}-non-class.csv
                else
                    ${dir}/build/drivers/osqp-bench -i $d -e $eps >> ${dir}/logs/perf_data/osqp-e${eps}-non-class.csv
                fi

                echo "Running OSQP-polished ..."
                if [ ! -f ${dir}/logs/perf_data/osqp-polished-e${eps}-non-class.csv ]; then
                    ${dir}/build/drivers/osqp-bench -i $d -e $eps -v polished -d 1 > ${dir}/logs/perf_data/osqp-polished-e${eps}-non-class.csv
                else
                    ${dir}/build/drivers/osqp-bench -i $d -e $eps -v polished >> ${dir}/logs/perf_data/osqp-polished-e${eps}-non-class.csv
                fi
            fi
        fi
    done
done

cd $dir

var="$(python -c 'import sys; print(sys.version_info[0])')"
if [[ $var == 2 ]]; then
    python3 all.py logs/perf_data
else
    python all.py logs/perf_data
fi