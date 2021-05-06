#!/bin/bash -l

# load appropriate modules
module load gcc/8.3.0
module load openmpi/2.1.2

# parameters
N="2000"
maxiter="1000"
chunksizes="1 3 10 32 100 316 1000 3162 10000 31623 100000 316228 1000000"

# mandelbrot basic
sbatch --export="N=${N},maxiter=${maxiter}" mandelbrot.slurm
echo

# mandelbrot static decomposition
sbatch --export="N=${N},maxiter=${maxiter}" mandelbrot-static.slurm
echo

# mandelbrot master-worker
for chunksize in ${chunksizes} ; do
    sbatch --export="N=${N},maxiter=${maxiter},chunksize=${chunksize}" \
           mandelbrot-master_worker.slurm
    echo
done

# mandelbrot cyclic decomposition
sbatch --export="N=${N},maxiter=${maxiter}" mandelbrot-cyclic.slurm
echo
