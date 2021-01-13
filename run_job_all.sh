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

# load all modules
module load cmake
module load gcc
module load intel
module load python
module load scipy-stack

./all.sh ../SMP_Repository $MKLROOT /home/zjming1/metis-5.1.0/
./all.sh ../SMP_Repository2 $MKLROOT /home/zjming1/metis-5.1.0/
./all_plot.sh