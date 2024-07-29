#!/bin/bash
#SBATCH -J m20U6.59
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --exclude=master,node[2,4,6-7,25-27]
#SBATCH -o log.out
#SBATCH -e log.err

program=./main.exe

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/opt/intel2019/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin/

mpiexec -np $SLURM_NPROCS -bootstrap ssh $program 
 
