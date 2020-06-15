#!/bin/bash


DATASET=SMP_Repository/
BUILDIR=build/

if [ "$#" -ge 2 ]; then
DATASET=$1
BUILDIR=$2
fi
eps=-3
echo "Running solvers in $BUILDIR for QP problems in $DATASET ..."


#for eps in {-3,-6,-9}; do
 echo "Running NASOQ-Fixed ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps> logs/nasoq-fixed-e${eps}.csv

 echo "Running NASOQ-Tuned ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v tuned"> logs/nasoq-tuned-e${eps}.csv

 echo "Running customized NASOQ ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET $eps "-v predet -r 0"> logs/nasoq-custom-e${eps}.csv

 echo "Running OSQP ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET $eps> logs/osqp-e${eps}.csv

 echo "Running OSQP-polished ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET $eps "-v polished"> logs/osqp-polished-e${eps}.csv

 echo "Running Gurobi ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/drivers/gurobi-bench $DATASET $eps > logs/gurobi-e${eps}.csv

 echo "Running Mosek ..."
 bash scripts/NASOQ_bench.sh $BUILDIR/drivers/mosek-bench $DATASET $eps > logs/mosek-e${eps}.csv
#done
sed -e "/Academic license - for non-commercial use only/d" -i logs/gurobi-*
echo "CSV files are generated in logs/"
echo "Plotting ..."
cd scripts/python_scripts/;
#for eps in {-3,-6,-9}; do
 python graph_generator.py -d ../../logs/ -s $eps
#done
rm -f *.txt
