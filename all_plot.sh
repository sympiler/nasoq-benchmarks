#!/bin/bash
export DISPLAY=unix$DISPLAY

var="$(python -c 'import sys; print(sys.version_info[0])')"
if [[ $var == 2 ]]; then
    python3 all.py logs/perf_data
else
    python all.py logs/perf_data
fi