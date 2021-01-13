#!/bin/bash
#SBATCH --mem=80GB

# load all modules
module load cmake
module load gcc
module load intel
module load python
module load scipy-stack

sbatch all.sh ../SMP_Repository $MKLROOT /home/zjming1/metis-5.1.0/
sbatch all.sh ../SMP_Repository2 $MKLROOT /home/zjming1/metis-5.1.0/
sbatch all_plot.sh