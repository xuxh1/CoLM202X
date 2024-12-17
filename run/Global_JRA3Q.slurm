#!/bin/bash
#SBATCH --job-name=JRA3Q_1129
#SBATCH -p first
#SBATCH -w tms04
#SBATCH -o JRA3Q_1129-o.out
#SBATCH -e JRA3Q_1129-e.out
#SBATCH -N 1 
#SBATCH -n 48
#SBATCH --ntasks-per-node=48
##SBATCH --exclusive
#SBATCH -t 2-24:00
##SBATCH --mem=1000G

module load intel-env

#which mpiifx
#which mpiifort
#which ifort

ulimit -s unlimited

echo "作业开始时间: $(date)"
echo "运行的节点: $(hostname)"

mpirun  -n 80  ./mksrfdata.x Global_JRA3Q.nml &> Global_JRA3Q_1129-o.log
mpirun  -n 80  ./mkinidata.x Global_JRA3Q.nml &>> Global_JRA3Q_1129-o.log
#mpirun  -np 1 ./colm.x Global_JRA3Q_cama.nml : -np 79 ./colm.x Global_JRA3Q_cama.nml &>> Global_cama_0.25_JRA3Q_1120-o.log
# /stu01/lianghb21/CoLM202X/2023student/lib/mpich/mpich-4.1.2/mpich-412/bin/mpirun  -np 1 ./colm.x Global_JRA3Q_cama.nml : -np 79 ./colm.x Global_JRA3Q_cama.nml &>> Global_cama_0.25_JRA3Q-o.log

echo "作业结束时间: $(date)"

2
