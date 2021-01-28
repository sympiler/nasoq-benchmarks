#!/bin/bash
# this command eneable plotting in the computing node
export DISPLAY=unix$DISPLAY

# call all.py to merge and generate plots
var="$(python -c 'import sys; print(sys.version_info[0])')"
if [[ $var == 2 ]]; then
    python3 all.py logs/perf_data
else
    python all.py logs/perf_data
fi