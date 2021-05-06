#!/bin/bash -l

# load appropriate modules
module load gcc/8.3.0
module load openmpi/2.1.2

# parameters
N="2000"
maxiter="1000"
chunksizes="10 100 1000 10000 100000 1000000"

# mandelbrot basic

# mandelbrot static decomposition

# mandelbrot master-worker
for chunksize in ${chunksizes} ; do
    sbatch --export="N=${N},maxiter=${maxiter},chunksize=${chunksize}" \
           mandelbrot_master_worker.slurm
done

# mandelbrot cyclic decomposition
