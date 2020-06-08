#!/bin/bash
DATASET=/home/kazem/SMP_Repository/


## NASOQ
#bash scripts/NASOQ_bench.sh nasoq/NASOQ-BIN $DATASET -3> logs/nasoq-fixed-e3.csv
bash scripts/NASOQ_bench.sh nasoq/NASOQ-BIN $DATASET -6> logs/nasoq-fixed-e6.csv
bash scripts/NASOQ_bench.sh nasoq/NASOQ-BIN $DATASET -9> logs/nasoq-fixed-e9.csv
