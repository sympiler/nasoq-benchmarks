#!/bin/bash
#SBATCH --mem=80GB

sbatch all.sh ../SMP_Repository $MKLROOT /home/zjming1/metis-5.1.0/
sbatch all.sh ../SMP_Repository2 $MKLROOT /home/zjming1/metis-5.1.0/