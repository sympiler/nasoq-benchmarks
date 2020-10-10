#!/bin/bash


DATASET=SMP_Repository/
BUILDIR=build/

if [ "$#" -ge 2 ]; then
DATASET=$1
BUILDIR=$2
fi

diag_perturb=-8
max_iter=3

mkdir -p ${DATASET}_test_settings_csvs

eps_lst='-3 -6'

for eps in $eps_lst
do
    echo "runing NASOQ-Fixed"
    if [ ! -f ${DATASET}_tune_csvs/nasoq-fixed-eps${eps}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb} -r ${max_iter}" \
            > ${DATASET}_test_settings_csvs/nasoq-fixed-eps${eps}.csv
    fi

    echo "runing NASOQ-Tuned"
    if [ ! -f ${DATASET}_tune_csvs/nasoq-tuned-eps${eps}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb} -r ${max_iter} -v tuned" \
            > ${DATASET}_test_settings_csvs/nasoq-tuned-eps${eps}.csv
    fi

    echo "runing NASOQ-Custom"
    if [ ! -f ${DATASET}_tune_csvs/nasoq-custom-eps${eps}.csv ]; then
            bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-p ${diag_perturb} -r ${max_iter} -v predet" \
            > ${DATASET}_test_settings_csvs/nasoq-custom-eps${eps}.csv
    fi
done