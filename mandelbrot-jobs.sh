#!/bin/bash -l

# load appropriate modules
module load gcc/8.3.0
module load openmpi/2.1.2

# parameters
N="2000"
maxiter="1000"
chunksizes="10 100 1000 10000 100000 1000000"

# compilation
gfortran src/mandelbrot.f90 -o bin/mandelbrot

mpifort src/mandelbrot_static.f90 -o bin/mandelbrot_static

mpifort src/mandelbrot_master_worker.f90 -o bin/mandelbrot_master_worker

# mandelbrot basic
sbatch --export="N=${N},maxiter=${maxiter}" mandelbrot.slurm

# mandelbrot static decomposition
sbatch --export="N=${N},maxiter=${maxiter}" mandelbrot_static.slurm

# mandelbrot master-worker
for chunksize in ${chunksizes} ; do
    sbatch --export="N=${N},maxiter=${maxiter},chunksize=${chunksize}" \
           mandelbrot_master_worker.slurm
done

# mandelbrot cyclic decomposition
