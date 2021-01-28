#!/bin/bash
#SBATCH --job-name="NASOQ_"
#SBATCH --output="NASOQ.%j.%N.out"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --export=ALL
#SBATCH -t 23:59:00
#SBATCH --mem=80GB
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
export OMP_NUM_THREADS=8

# this script invoke tests in the repo

# load all modules
module load cmake
module load gcc
module load intel
module load python
module load scipy-stack

# invoke all testing process
# please modify here based on where you installed MKL and metis
# generate performance data at first, and then make tests and generate plots

# in the format:
# ./all.sh <dataset> <MKL path> <METIS path>
# ...
# ./all_plot.sh: make plots
./all.sh ../SMP_Repository $MKLROOT /home/zjming1/metis-5.1.0/
./all.sh ../SMP_Repository2 $MKLROOT /home/zjming1/metis-5.1.0/
./all_plot.sh