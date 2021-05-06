#!/bin/bash -l

# SLURM details
#SBATCH --account=courses0100
#SBATCH --reservation=courses0100
#SBATCH --job-name=mandelbrot
#SBATCH --time=00:10:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --tasks-per-node=10
#SBATCH --cpus-per-task=1

# load appropriate modules
module load gcc/8.3.0

# parameters
# N="2000"
# maxiter="1000"

# compilation
gfortran src/mandelbrot.f90 -o bin/mandelbrot

# mandelbrot basic
srun -n 1 bin/mandelbrot ${N} ${maxiter}
