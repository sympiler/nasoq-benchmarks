# nasoq_benchmarks


## Building benchmark
OSQP and NASOQ are already in the repository however, 
you need to get the license for other tools, i.e., Gurobi, 
Mosek, and QL. 

### Requirements
- NASOQ needs MKL Pardiso and METIS.
- OSQP needs MKL library for fast solve time.
- MOSEK and Gurobi are commercial solver.

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DMETI ..
make
```

## Running the benchmark
You first need to download the dataset and then run the all compiled 
solvers.

### Downloading the dataset
```bash
wget 
```
### Running solvers
A script is provided that run all tools for specified dataset.
To use the script emit the following command:
```bash
bash run_all path/to/dataset/ path/to/build/directory/
```

The script finds all sparse QP problems stored in `yml`  and 
run all solvers for each of them for three different ranges of 
accuracy, i.e., 1e-3, 1e-6, and 1e-9. 

