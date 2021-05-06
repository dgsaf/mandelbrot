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
# srun -n 1 bin/mandelbrot ${N} ${maxiter}

# mandelbrot static decomposition
# srun -n 10 bin/mandelbrot_static ${N} ${maxiter}

# mandelbrot master-worker
for chunksize in ${chunksizes} ; do
    srun -n 10 bin/mandelbrot_master_worker ${N} ${maxiter} ${chunksize}
done

# mandelbrot cyclic decomposition
