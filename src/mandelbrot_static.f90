program mandelbrot_static

  use mpi

  complex :: z, kappa
  integer :: green, blue, i, j, k, loop

  integer , parameter :: N = 2000
  integer , parameter :: maxiter = 1000
  real :: x(0:N*N-1)

  ! MPI variables.
  integer :: proc_id, n_proc, err
  integer :: chunksize
  integer , allocatable :: loop_min(:), loop_max(:), chunksize_proc(:)
  real , allocatable :: x_proc(:)
  integer :: p

  ! MPI timing variables.
  double precision :: times(1:7)
  double precision :: time_setup, time_comp, time_wait, time_comm, time_total

  ! MPI initialisation.
  call MPI_INIT(err)

  times(1) = MPI_WTIME()

  call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

  ! Determine (default) chunk size.
  chunksize = ceiling(real(N*N) / real(n_proc))

  ! Determine loop bounds and chunk size for each process.
  allocate(loop_min(0:n_proc-1))
  allocate(loop_max(0:n_proc-1))
  allocate(chunksize_proc(0:n_proc-1))

  do p = 0, n_proc - 1
    loop_min(p) = max(0, chunksize * p)
    loop_max(p) = min(N*N-1, chunksize * (p + 1) - 1)
    chunksize_proc(p) = loop_max(p) - loop_min(p) + 1
  end do

  ! Allocate work array for given process.
  allocate(x_proc(0:chunksize_proc(proc_id)-1))

  times(2) = MPI_WTIME()

  ! Mandelbrot calculation.
  ! Modified to loop only over a certain subset of indexes (due to utilising
  ! static decomposition).
  do loop = loop_min(proc_id), loop_max(proc_id)
    ! i varies from 0 to N-1
    i = int(loop/N)

    ! j varies from 0 to N-1
    j = mod(loop, N)
    kappa = cmplx((4.0*(i-N/2)/N), (4.0*(j-N/2)/N))

    k = 1
    z = kappa
    do while (abs(z) <= 2 .and. k < maxiter)
      k = k+1
      z = z*z + kappa
    end do

    ! x_proc is indexed to ensure all x_proc(0:chunksize_proc(proc_id)) are
    ! determined.
    x_proc(loop-loop_min(proc_id)) = log(real(k))/log(real(maxiter))
  end do

  times(3) = MPI_WTIME()

  ! MPI wait for all process to finish
  call MPI_BARRIER(MPI_COMM_WORLD, err)

  times(4) = MPI_WTIME()

  ! MPI gather the work array from each process to the root process.
  call MPI_GATHERV(x_proc, chunksize_proc(proc_id), MPI_REAL, x, &
      chunksize_proc, loop_min, MPI_REAL, 0, MPI_COMM_WORLD, err)

  times(5) = MPI_WTIME()

  ! Timing analysis.
  time_setup = times(2) - times(1)
  time_comp = times(3) - times(2)
  time_wait = times(4) - times(3)
  time_comm = times(5) - times(4)
  time_total = times(5) - times(1)

  write (*, '(a, i3, a, f8.5, a, f8.5, a, f8.5, a, f8.5, a)') &
      "timing for process: ", proc_id, " (", &
      100d0*time_setup/time_total, " %, ", &
      100d0*time_comp/time_total, " %, ", &
      100d0*time_wait/time_total, " %, ", &
      100d0*time_comm/time_total, " %)"

  ! write (*, *) &
  !     "timing for process: ", proc_id, NEW_LINE('a'), &
  !     "  setup:         ", time_setup, NEW_LINE('a'), &
  !     "  computation:   ", time_comp, NEW_LINE('a'), &
  !     "  waiting:       ", time_wait, NEW_LINE('a'), &
  !     "  communicating: ", time_comm, NEW_LINE('a'), &
  !     "  total: ", time_total

  ! Writing data to file (only for root process).
  if (proc_id == 0) then

    write(*, *) "Writing mandelbrot_static.ppm"

    open(7, file="mandelbrot_static.ppm", status="unknown")

    write(7, 100) "P3", N, N, 255

    do loop = 0, N*N-1
      if (x(loop) < 0.5) then
        green = 2.0*x(loop)*255
        write(7, 110) 255-green, green, 0
      else
        blue = 2.0*x(loop)*255 - 255
        write(7, 110) 0, 255-blue, blue
      end if
    end do

100 format(A2, /, I4, I5, /, I3)
110 format(I3, /, I3, /, I3)

  end if

  ! mpi finalisation
  call MPI_FINALIZE(err)

  close(7)
end program mandelbrot_static