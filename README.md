# nasoq_benchmarks


## Building benchmark
OSQP and NASOQ are already in the repository however, 
you need to get the license for other tools, i.e., Gurobi, 
Mosek, and QL. 

### Requirements
- NASOQ needs MKL Pardiso and METIS. CMake 
downloads the code from NASOQ repository.
- OSQP needs MKL library for fast solve time. CMake downloads 
the source from OSQP repository. 
- MOSEK : Put the MOSEK license in 
```bash
export MOSEKROOT=path/to/mosek/8/tools/platform/linux64x86/
```
- Gurobi: Put the Gurobi license in the home directory and 
export the path that gurobi is installed, for example:
```bash
export GUROBIROOT=/opt/gurobi811/linux64/
```

### Build
First clone the repo and do not forget ```-recursive``` option:
```bash
git clone --recursive https://github.com/sympiler/nasoq-benchmarks
```

Before building the benchmark, you need to know the path 
to METIS and MKL library (CMake can find MKL if it is 
installed in the default location or is in the path). For MKL 
you need to set `MKLROOT` and for METIS you need to set 
`METIS_ROOT_PATH`. Following shows how these variables are set using 
cmake flags.

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DMETIS_ROOT_PATH=path/to/metis-5.1.0/build/Linux-x86_64/ ..
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

