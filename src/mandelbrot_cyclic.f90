!> mandelbrot_cyclic
!
!  Adaption of mandelbrot program to use MPI with a cyclic decomposition of the
!  data.
!
!  For a NxN grid (represented in a 1D array; that is, {0, .., N*N - 1}),
!  with n_proc processes.
!
!  The default chunksize, the size of the data subset to be assigned to each
!  process, is defined to be the smallest integer such that
!  > chunksize * n_proc >= N*N.
!  This is so that even in the case where n_proc doesn't divide N*N, every
!  data point is guaranteed to be covered by a process. However, we note that
!  the remainder of this program has not been adapted to work in cases where
!  n_proc doesn't divide N*N. May lead to segfault errors in such a case.
!
!  Each process, p = 0, .., n_proc-1, is assigned a section of the data
!  {p + (n_proc * k) , k = 0, ..., chunksize-1 }, to perform the mandelbrot
!  calculation on.
!
!  The results of each process are then transposed to recover the
!  block-sequential, rather than cyclic, data, which is then gathered at the
!  root.
program mandelbrot_cyclic

  use mpi

  complex :: z, kappa
  integer :: green, blue, i, j, k, loop

  integer :: N, maxiter
  real , allocatable :: x(:)

  ! MPI variables.
  !  proc_id         ID of current MPI process.
  !  n_proc          Number of MPI processes.
  !  err             Stores error code for MPI calls.
  !  chunksize       The default chunksize; the maximum possible number of loop
  !                  iterations assigned to each process. Defined to be the
  !                  smallest integer such that
  !                  that chunksize * n_proc >= N*N.
  !                  Note that this program has not been adapted to work for
  !                  cases where n_proc does not divide N*N.
  !  x_proc          A work array, local to each process, which will yield the
  !                  mandelbrot data for the data subset assigned to this
  !                  process. Will have size chunksize.
  !  x_cyclic        A work array, for the root process, which will gather the
  !                  cyclic data, and from which x(:) can be recovered
  !  proc            A counter variable for looping over processes.
  !  l               A counter variable used for clarity when performing cyclic
  !                  loops.
  integer :: proc_id, n_proc, err
  integer :: chunksize
  real , allocatable :: x_proc(:), x_cyclic(:)
  integer :: proc, l

  ! Timing variables.
  !  times       An array storing time markers used to determine the following
  !              timing variables.
  !  time_setup  Time taken for this process to setup MPI variables for
  !              partitioning data.
  !  time_comp   Time taken for this process to perform mandelbrot calculations
  !              for its given data subset.
  !  time_wait   Time this process spends waiting while other processes finish
  !              performing their calculations.
  !  time_comm   Time taken for this process to communicate its data subset
  !              to the root process.
  !  time_total  Time taken overall.
  double precision :: times(1:5)
  double precision :: time_setup, time_comp, time_wait, time_comm, time_total

  ! Read command line arguments.
  call read_input(N, maxiter)

  allocate(x(0:N*N-1))

  ! MPI initialisation.
  call MPI_INIT(err)

  times(1) = MPI_WTIME()

  call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

  ! Determine (default) chunk size.
  chunksize = ceiling(real(N*N) / real(n_proc))

  ! Allocate work array for given process.
  allocate(x_proc(0:chunksize-1))
  allocate(x_cyclic(0:N*N-1))

  times(2) = MPI_WTIME()

  ! Mandelbrot calculation.
  ! Modified to loop only over a certain subset of indexes (due to utilising
  ! cyclic decomposition).
  do l = 0, chunksize-1
    loop = (n_proc * l) + proc_id

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

    x_proc(l) = log(real(k))/log(real(maxiter))
  end do

  times(3) = MPI_WTIME()

  ! MPI wait for all process to finish.
  call MPI_BARRIER(MPI_COMM_WORLD, err)

  times(4) = MPI_WTIME()

  ! MPI gather the work array from each process to the root process.
  call MPI_GATHER(x_proc, chunksize, MPI_REAL, x_cyclic, chunksize, MPI_REAL, &
      0, MPI_COMM_WORLD, err)

  ! Extract the original data from the cyclically distributed data
  if (proc_id == 0) then

    do loop = 0, N*N-1
      proc = mod(loop, n_proc)
      l = (loop - proc) / n_proc

      x(loop) = x_cyclic((chunksize*proc) + l)
    end do

  end if

  times(5) = MPI_WTIME()

  ! Timing analysis.
  time_setup = times(2) - times(1)
  time_comp = times(3) - times(2)
  time_wait = times(4) - times(3)
  time_comm = times(5) - times(4)
  time_total = times(5) - times(1)

  if (proc_id == 0) then
    write (*, *) &
        "timing for MPI cyclic code:", NEW_LINE('a'), &
        "  total: ", time_total

    write (*, *) "time spent working/waiting/communicating:"
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, err)

  write (*, '(a, i1, a, f5.2, a, f5.2, a, f5.2, a)') &
      "  ", proc_id, ": ", &
      100.0*time_comp/time_total, " % / ", &
      100.0*time_wait/time_total, " % / ", &
      100.0*time_comm/time_total, " %"

  ! Write timing data to an output file.
  call write_timing_data (N, maxiter, n_proc, proc_id, &
      time_setup, time_comp, time_wait, time_comm, time_total)

  call MPI_BARRIER(MPI_COMM_WORLD, err)

!   ! Writing data to file (only done by root process).
!   if (proc_id == 0) then

!     write (*, *) "Writing mandelbrot_cyclic.ppm"

!     open(7, file="mandelbrot_cyclic.ppm", status="unknown")

!     write(7, 100) "P3", N, N, 255

!     do loop = 0, N*N-1
!       if (x(loop) < 0.5) then
!         green = 2.0*x(loop)*255
!         write(7, 110) 255-green, green, 0
!       else
!         blue = 2.0*x(loop)*255 - 255
!         write(7, 110) 0, 255-blue, blue
!       end if
!     end do

! 100 format(A2, /, I4, I5, /, I3)
! 110 format(I3, /, I3, /, I3)

!   end if

  ! MPI finalisation.
  call MPI_FINALIZE(err)

  close(7)

contains

  ! Read in the value of N, maxiter from command line arguments
  subroutine read_input (N, maxiter)
    integer , intent(out) :: N, maxiter
    integer :: num_args
    character(len=20) :: arg

    num_args = command_argument_count()

    if (num_args < 2) then
      write (*, *) "Usage: <N> <maxiter>"
    end if

    if (num_args >= 1) then
      call get_command_argument(1, arg)
      read (arg, *) N
    else
      write (*, *) "<N> not specified. Using default value, N = 2000."
      N = 2000
    end if

    if (num_args >= 2) then
      call get_command_argument(2, arg)
      read (arg, *) maxiter
    else
      write (*, *) &
          "<maxiter> not specified. Using default value, maxiter = 1000."
      maxiter = 1000
    end if

  end subroutine read_input

  ! Write the timing data (for the given parameters: N, maxiter, chunksize,
  ! n_proc), for a given process.
  !
  ! The timing data includes the time taken spent: setting up, communicating,
  ! performing computations, waiting, and the total time spent.
  !
  ! The filename is defined by (N, maxiter, n_proc, proc_id)
  subroutine write_timing_data (N, maxiter, n_proc, proc_id, &
      time_setup, time_comp, time_wait, time_comm, time_total)
    integer , intent(in) :: N, maxiter, proc_id
    double precision , intent(in) :: time_setup, time_comp, time_wait, &
        time_comm, time_total
    character(len=1000) :: timing_file
    character(len=20) :: str_N, str_maxiter, str_n_proc, str_proc_id
    integer :: file_unit

    ! Construct timing filename to be of the form:
    ! "output/timing.cyclic.N=<N>.maxiter=<maxiter>.n_proc=<n_proc>\
    ! .proc_id=<proc_id>.dat"
    write (str_N, *) N
    write (str_maxiter, *) maxiter
    write (str_n_proc, *) n_proc
    write (str_proc_id, *) proc_id

    write (timing_file, *) &
        "output/timing.cyclic.", &
        "N-", trim(adjustl(str_N)), ".", &
        "maxiter-", trim(adjustl(str_maxiter)), ".", &
        "n_proc-", trim(adjustl(str_n_proc)), ".", &
        "proc_id-", trim(adjustl(str_proc_id)), ".dat"

    ! Append the timing data to the data file
    file_unit = 10 + proc_id

    open (file_unit, file=trim(adjustl(timing_file)), action="write")

    write (file_unit, *) time_setup, time_comp, time_wait, &
        time_comm, time_total

    close (file_unit)

  end subroutine write_timing_data

end program mandelbrot_cyclic
