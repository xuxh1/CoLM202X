#!/bin/bash
#SBATCH -J US_cama
#SBATCH -p second
#SBATCH -w tms04
#SBATCH -N 1 
#SBATCH -n 48
#SBATCH --ntasks-per-node=80
##SBATCH --mem=1000G
#SBATCH -o US_cama-e.%j
#SBATCH -e US_cama-o.%j
##SBATCH --exclusive
#SBATCH -t 5-24:00

ulimit -s unlimited

# mpirun -n 48 ./mksrfdata.x test.nml
# mpirun -n 48 ./mkinidata.x test.nml
# mpirun -n 48 ./colm.x test.nml

mpirun -n 80 ./mksrfdata.x US_JRA55_cama.nml
mpirun -n 80 ./mkinidata.x US_JRA55_cama.nml
mpirun -n 1 ./colm.x US_JRA55_cama.nml : -n 79 ./colm.x US_JRA55_cama.nml

# mpirun -n 12 ./mksrfdata.x test_wa.nml
# mpirun -n 12 ./mkinidata.x test_wa.nml
# mpirun -n 12 ./colm.x test_wa.nml

