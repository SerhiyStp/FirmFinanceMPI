#!/bin/bash
#SBATCH -n 100                # requested nodes
#SBATCH -p parallel-short     # requested queue
#SBATCH -t 300                # maximum runtime in minutes

module load mvapich2/2.0-gcc-4.9.1
#module load mvapich2/2.1-gcc-5.1.0
#module load mvapich2/2.2b-pgi-15.4

srun -n 100 ./firm_finance.exe
