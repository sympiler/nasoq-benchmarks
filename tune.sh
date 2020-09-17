#!/bin/bash


DATASET=SMP_Repository/
BUILDIR=build/

if [ "$#" -ge 2 ]; then
DATASET=$1
BUILDIR=$2
fi

eps_lst='-3 -6'
max_iter_lst='0 1 2 3 4 5 10 20'
stop_tol_lst='-13 -15 -16 -17'
diag_perturb_lst='-6 -7 -8 -9 -10 -11 -12'

default_max_iter=0
default_stop_tol=-13
default_diag_perturb=-6

mkdir -p tune_csvs
mkdir -p tune_csvs/max_iter
mkdir -p tune_csvs/stop_tol
mkdir -p tune_csvs/diag_perturb

for eps in $eps_lst do
    for max_iter in $max_iter_lst do
        if [ ! -f tune_csvs/max_iter/nasoq-fixed-eps${eps}_max_iter${max_iter}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${default_diag_perturb} -r ${max_iter} -t ${default_stop_tol}" \
            > tune_csvs/max_iter/nasoq-fixed-eps${eps}_max_iter${max_iter}.csv
        fi

        if [ ! -f tune_csvs/max_iter/nasoq-tuned-eps${eps}_max_iter${max_iter}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${default_diag_perturb} -r ${max_iter} -t ${default_stop_tol} -v tuned" \
            > tune_csvs/max_iter/nasoq-tuned-eps${eps}_max_iter${max_iter}.csv
        fi

        if [ ! -f tune_csvs/max_iter/nasoq-custom-eps${eps}_max_iter${max_iter}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${default_diag_perturb} -r ${max_iter} -t ${default_stop_tol} -v predet" \
            > tune_csvs/max_iter/nasoq-custom-eps${eps}_max_iter${max_iter}.csv
        fi
    done

    for stop_tol in $stop_tol_lst do
        if [ ! -f tune_csvs/stop_tol/nasoq-fixed-eps${eps}_stop_tol${stop_tol}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${default_diag_perturb} -r ${default_max_iter} -t ${stop_tol}" \
            > tune_csvs/stop_tol/nasoq-fixed-eps${eps}_stop_tol${stop_tol}.csv
        fi

        if [ ! -f tune_csvs/stop_tol/nasoq-tuned-eps${eps}_stop_tol${stop_tol}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${default_diag_perturb} -r ${default_max_iter} -t ${stop_tol} -v tuned" \
            > tune_csvs/stop_tol/nasoq-tuned-eps${eps}_stop_tol${stop_tol}.csv
        fi

        if [ ! -f tune_csvs/stop_tol/nasoq-custom-eps${eps}_stop_tol${stop_tol}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${default_diag_perturb} -r ${default_max_iter} -t ${stop_tol} -v predet" \
            > tune_csvs/stop_tol/nasoq-custom-eps${eps}_stop_tol${stop_tol}.csv
        fi
    done

    for diag_perturb in $diag_perturb_lst do
        if [ ! -f tune_csvs/diag_perturb/nasoq-fixed-eps${eps}_diag_perturb${diag_perturb}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb} -r ${default_max_iter} -t ${default_stop_tol}" \
            > tune_csvs/diag_perturb/nasoq-fixed-eps${eps}_diag_perturb${diag_perturb}.csv
        fi

        if [ ! -f tune_csvs/diag_perturb/nasoq-tuned-eps${eps}_diag_perturb${diag_perturb}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb} -r ${default_max_iter} -t ${default_stop_tol} -v tuned" \
            > tune_csvs/diag_perturb/nasoq-tuned-eps${eps}_diag_perturb${diag_perturb}.csv
        fi

        if [ ! -f tune_csvs/diag_perturb/nasoq-custom-eps${eps}_diag_perturb${diag_perturb}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb} -r ${default_max_iter} -t ${default_stop_tol} -v predet" \
            > tune_csvs/diag_perturb/nasoq-custom-eps${eps}_diag_perturb${diag_perturb}.csv
        fi
    done
done