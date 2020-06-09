#!/bin/bash
DATASET=/home/kazem/SMP_Repository/alligator/
BUILDIR=/home/kazem/development/nasoq-benchmarks/cmake-build-debug/


## NASOQ
bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET -3> logs/nasoq-fixed-e3.csv
bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET -6> logs/nasoq-fixed-e6.csv
bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET -9> logs/nasoq-fixed-e9.csv

## NASOQ-Tuned
bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET -3 "-v tuned"> logs/nasoq-tuned-e3.csv
bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET -6 "-v tuned"> logs/nasoq-tuned-e6.csv
bash scripts/NASOQ_bench.sh $BUILDIR/nasoq/NASOQ-BIN $DATASET -9 "-v tuned"> logs/nasoq-tuned-e9.csv

## OSQP
bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET -3> logs/osqp-e3.csv
bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET -6> logs/osqp-e6.csv
bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET -9> logs/osqp-e9.csv

## OSQP-Polished
bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET -3 "-v polished"> logs/osqp-polished-e3.csv
bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET -6 "-v polished"> logs/osqp-polished-e6.csv
bash scripts/NASOQ_bench.sh $BUILDIR/drivers/osqp-bench $DATASET -9 "-v polished"> logs/osqp-polished-e9.csv
