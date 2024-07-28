#!/bin/bash
#SBATCH -J Validate
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o log.out
#SBATCH -e log.err

program=./main.exe

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/opt/intel2019/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin/

mpiexec -np $SLURM_NPROCS -bootstrap ssh $program 
 
