#!/bin/bash -l

# load appropriate modules
module load gcc/8.3.0
module load openmpi/2.1.2

# parameters
N="2000"
maxiter="1000"
chunksizes="1 3 10 32 100 316 1000 3162 10000 31623 100000 316228 1000000"

# compilation
gfortran src/mandelbrot.f90 -o bin/mandelbrot

mpifort src/mandelbrot_static.f90 -o bin/mandelbrot_static

mpifort src/mandelbrot_master_worker.f90 -o bin/mandelbrot_master_worker

mpifort src/mandelbrot_cyclic.f90 -o bin/mandelbrot_cyclic

# mandelbrot basic
srun -n 1 bin/mandelbrot ${N} ${maxiter}
echo

# mandelbrot static decomposition
srun -n 10 bin/mandelbrot_static ${N} ${maxiter}
echo

# mandelbrot master-worker
for chunksize in ${chunksizes} ; do
    srun -n 10 bin/mandelbrot_master_worker ${N} ${maxiter} ${chunksize}
    echo
done

# mandelbrot cyclic decomposition
srun -n 10 bin/mandelbrot_cyclic ${N} ${maxiter}
echo
