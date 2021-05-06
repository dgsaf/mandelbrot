#!/bin/bash -l

# SLURM details
#SBATCH --account=courses0100
#SBATCH --reservation=courses0100
#SBATCH --job-name=mandelbrot-static
#SBATCH --time=00:10:00
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --tasks-per-node=10
#SBATCH --cpus-per-task=1

# load appropriate modules
module load gcc/8.3.0
module load openmpi/2.1.2

# parameters
# N="2000"
# maxiter="1000"

# compilation
mpifort src/mandelbrot_static.f90 -o bin/mandelbrot_static

# mandelbrot static decomposition
srun -n 10 bin/mandelbrot_static ${N} ${maxiter}
