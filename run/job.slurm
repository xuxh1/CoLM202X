#!/bin/bash
#SBATCH -J Siberia_CH4_1
#SBATCH -p first
#SBATCH -w tms03
#SBATCH -N 1 
#SBATCH -n 48
#SBATCH --ntasks-per-node=48
##SBATCH --mem=1000G
#SBATCH -o test_wa.%j
#SBATCH -e test_wa.%j
##SBATCH --exclusive
#SBATCH -t 2-24:00

ulimit -s unlimited

# mpirun -n 48 ./mksrfdata.x test.nml
# mpirun -n 48 ./mkinidata.x test.nml
# mpirun -n 48 ./colm.x test.nml

mpirun -n 48 ./mksrfdata.x test.nml
mpirun -n 48 ./mkinidata.x test.nml
mpirun -n 48 ./colm.x test.nml

# mpirun -n 12 ./mksrfdata.x test_wa.nml
# mpirun -n 12 ./mkinidata.x test_wa.nml
# mpirun -n 12 ./colm.x test_wa.nml

