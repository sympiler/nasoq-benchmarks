#!/bin/bash
#SBATCH --mem=80GB

sbatch all.sh ../SMP_Repository $MKLROOT /home/zjming1/metis-5.1.0/
sbatch all.sh ../SMP_Repository2 $MKLROOT /home/zjming1/metis-5.1.0/

var="$(python -c 'import sys; print(sys.version_info[0])')"
if [[ $var == 2 ]]; then
    python3 all.py logs/perf_data
else
    python all.py logs/perf_data
fi