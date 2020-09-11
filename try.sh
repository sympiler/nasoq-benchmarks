#!/bin/bash


DATASET=SMP_Repository/
BUILDIR=build/

if [ "$#" -ge 3 ]; then
DATASET=$1
BUILDIR=$2
eps=$3
fi
# eps=-3
echo "Running solvers in $BUILDIR for QP problems in $DATASET ..."


# for eps in {-3,-6}; do

    for diag_perturb in {-6, -9, -12}; do
        echo "Running NASOQ-Fixed ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb}"> logs/nasoq-fixed-e${eps}.csv

        echo "Running NASOQ-Tuned ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v tuned -p ${diag_perturb}"> logs/nasoq-tuned-e${eps}.csv

        echo "Running customized NASOQ ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v predet -r 0 -p ${diag_perturb}"> logs/nasoq-custom-e${eps}.csv

        cd scripts/python_scripts/;
        python graph_generator.py -d ../../logs/ -s $eps
        rm -f *.txt
        cd ../..
    done

    for stop_tol in {-13, -15, -17}; do
        echo "Running NASOQ-Fixed ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-t ${stop_tol}"> logs/nasoq-fixed-e${eps}.csv

        echo "Running NASOQ-Tuned ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v tuned -t ${stop_tol}"> logs/nasoq-tuned-e${eps}.csv

        echo "Running customized NASOQ ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v predet -r 0 -t ${stol_tol}"> logs/nasoq-custom-e${eps}.csv

        cd scripts/python_scripts/;
        python graph_generator.py -d ../../logs/ -s $eps
        rm -f *.txt
        cd ../..
    done

    for max_iter in {0, 5, 10}; do
        echo "Running NASOQ-Fixed ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-r ${max_iter}"> logs/nasoq-fixed-e${eps}.csv

        echo "Running NASOQ-Tuned ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v tuned -r ${max_iter}"> logs/nasoq-tuned-e${eps}.csv

        echo "Running customized NASOQ ..."
        bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v predet -r 0 -r ${max_iter}"> logs/nasoq-custom-e${eps}.csv

        cd scripts/python_scripts/;
        python graph_generator.py -d ../../logs/ -s $eps
        rm -f *.txt
        cd ../..
    done

# done

echo done
