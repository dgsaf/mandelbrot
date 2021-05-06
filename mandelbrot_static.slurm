#!/bin/bash -l

# SLURM details
#SBATCH --account=courses0100
#SBATCH --reservation=courses0100
#SBATCH --job-name=mandelbrot
#SBATCH --time=00:10:00
#SBATCH --export=NONE

# load appropriate modules
module load gcc/8.3.0
module load openmpi/2.1.2

# program execution
srun -n 10 mandelbrot_static ${N} ${maxiter}
