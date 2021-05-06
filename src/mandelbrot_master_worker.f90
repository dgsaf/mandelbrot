!> mandelbrot_master_worker
!!
!!  Adaption of mandelbrot program to use MPI with master_worker parallelism.
!!
!!  For a NxN grid (represented in a 1D array; that is, {0, .., N*N - 1}),
!!  with n_proc processes.
!!
program mandelbrot_master_worker

  use mpi

  complex :: z, kappa
  integer :: green, blue, i, j, k, loop

  integer :: N, maxiter
  real , allocatable :: x(:)

  ! MPI variables.
  !!  proc_id         ID of current MPI process.
  !!  n_proc          Number of MPI processes.
  !!  err             Stores error code for MPI calls.
  !!  chunksize       The number of loop iterations assigned to a process for
  !!                  one task.
  !!  n_tasks         The number of tasks.
  !!  loop_min        An array of the lower loop-iteration bounds for each
  !!                  task.
  !!  loop_max        An array of the upper loop-iteration bounds for each
  !!                  task.
  !!  x_proc          A work array, local to each process, which will yield the
  !!                  mandelbrot data for the data subset assigned to this
  !!                  process. Will have size chunksize_proc(proc_id).
  !!  task            A counter variable for looping over tasks.
  !!  proc            A counter variable for looping over processes.
  integer , parameter :: master_id = 0
  integer :: proc_id, n_proc, err, tag, request
  integer :: status(MPI_STATUS_SIZE)
  integer :: chunksize, n_tasks
  integer , allocatable :: loop_min(:), loop_max(:)
  real , allocatable :: x_task(:)
  integer :: proc, task
  integer :: proc_recv, task_recv
  logical :: all_tasks_distributed
  integer , allocatable :: task_ledger(:)
  integer , parameter :: no_task = 0

  ! Timing variables.
  !!  times       An array storing time markers used to determine the following
  !!              timing variables.
  !!  time_setup  Time taken for this process to setup MPI variables for
  !!              partitioning data.
  !!  time_comp   Time taken for this process to perform mandelbrot calculations
  !!              for its given data subset.
  !!  time_wait   Time this process spends waiting while other processes finish
  !               performing their calculations.
  !!  time_comm   Time taken for this process to communicate its data subset
  !!              to the root process.
  !!  time_total  Time taken overall.
  double precision :: times(1:7)
  double precision :: time_setup, time_comp, time_wait, time_comm, time_total

  ! Read command line arguments.
  call read_input(N, maxiter, chunksize)

  allocate(x(0:N*N-1))

  ! Initialise timing variables
  time_setup = 0d0
  time_comp = 0d0
  time_wait = 0d0
  time_comm = 0d0
  time_total = 0d0

  ! MPI initialisation.
  call MPI_INIT(err)

  times(1) = MPI_WTIME()

  call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

  ! Arbitrarily set tag = 0.
  tag = 0

  ! Determine number of tasks needed for given chunksize.
  n_tasks = ceiling(real(N*N) / real(chunksize))

  allocate(loop_min(1:n_tasks))
  allocate(loop_max(1:n_tasks))

  do task = 1, n_tasks
    loop_min(task) = max(0, chunksize * (task - 1))
    loop_max(task) = min(N*N-1, chunksize * task - 1)

    ! !debugging
    ! if (proc_id == master_id) then
    !   write (*, *) "task", task, "has bounds (", &
    !       loop_min(task), " , ", loop_max(task), ")"
    ! end if
  end do

  ! Default value for flag indicating all tasks have been handed out
  all_tasks_distributed = .false.

  ! Allocate work arrays.
  allocate(x_task(0:chunksize-1))

  ! Master-Worker Scheme
  if (proc_id == master_id) then
    ! Master process.

    ! Allocate process ledger (records which task a given process is working
    ! on).
    allocate(task_ledger(1:n_proc-1))
    task_ledger(:) = no_task

    times(2) = MPI_WTIME()
    time_setup = times(2) - times(1)

    ! Initialise task counter.
    task = 1

    ! Distribute initial tasks.
    do proc = 1, n_proc - 1

      call MPI_ISEND(all_tasks_distributed, 1, MPI_LOGICAL, proc, tag, &
          MPI_COMM_WORLD, request, err)

      if (all_tasks_distributed) then
        ! !debugging
        ! write (*, '(a, i3, a, i10)') &
        !     "master -> ", proc_recv, " : no more tasks"

        task_ledger(proc) = no_task
      else
        call MPI_ISEND(task, 1, MPI_INTEGER, proc, tag, MPI_COMM_WORLD, &
            request, err)

        ! !debugging
        ! write (*, '(a, i3, a, i10)') &
        !     "master -> ", proc, " : task", task

        task_ledger(proc) = task
        task = task + 1
      end if

      all_tasks_distributed = (task > n_tasks)

      ! !debugging
      ! write (*, *) task_ledger(:)
    end do

    times(3) = MPI_WTIME()
    time_comm = time_comm + times(3) - times(2)

    ! Loop until all tasks distributed.
    do while (.not. all_tasks_distributed)
      times(4) = MPI_WTIME()

      ! Receive completed task from worker.
      call MPI_RECV(x_task, chunksize, MPI_REAL, MPI_ANY_SOURCE, tag, &
          MPI_COMM_WORLD, status, err)

      proc_recv = status(MPI_SOURCE)
      task_recv = task_ledger(proc_recv)

      ! !debugging
      ! write (*, '(a, i3, a, i10)') &
      !     "master <- ", proc_recv, " : task", task_recv

      times(5) = MPI_WTIME()

      x(loop_min(task_recv):loop_max(task_recv)) = x_task(:)

      times(6) = MPI_WTIME()

      ! Tell worker there is more work to do.
      call MPI_ISEND(all_tasks_distributed, 1, MPI_LOGICAL, proc_recv, tag, &
          MPI_COMM_WORLD, request, err)

      ! Distribute new task.
      call MPI_ISEND(task, 1, MPI_INTEGER, proc_recv, tag, MPI_COMM_WORLD, &
          request, err)

      ! !debugging
      ! write (*, '(a, i3, a, i10)') &
      !     "master -> ", proc_recv, " : task", task

      task_ledger(proc_recv) = task
      task = task + 1

      all_tasks_distributed = (task > n_tasks)

      ! !debugging
      ! write (*, *) task_ledger(:)

      times(7) = MPI_WTIME()

      time_comp = time_comp + times(6) - times(5)
      time_comm = time_comm + (times(5) - times(4)) + (times(7) - times(6))
    end do

    ! Collect outstanding tasks.
    do while (any(task_ledger /= no_task))
      times(4) = MPI_WTIME()

      ! Receive completed task from worker.
      call MPI_RECV(x_task, chunksize, MPI_REAL, MPI_ANY_SOURCE, tag, &
          MPI_COMM_WORLD, status, err)

      proc_recv = status(MPI_SOURCE)
      task_recv = task_ledger(proc_recv)

      ! !debugging
      ! write (*, '(a, i3, a, i10)') &
      !     "master <- ", proc_recv, " : task", task_recv

      times(5) = MPI_WTIME()

      x(loop_min(task_recv):loop_max(task_recv)) = x_task(:)

      times(6) = MPI_WTIME()

      ! Tell worker there is no more work to do.
      call MPI_ISEND(all_tasks_distributed, 1, MPI_LOGICAL, proc_recv, tag, &
          MPI_COMM_WORLD, request, err)

      ! !debugging
      ! write (*, '(a, i3, a, i10)') &
      !     "master -> ", proc_recv, " : no more tasks"

      ! No more work to distribute.
      task_ledger(proc_recv) = no_task

      ! !debugging
      ! write (*, *) task_ledger(:)

      times(7) = MPI_WTIME()

      time_comp = time_comp + times(6) - times(5)
      time_comm = time_comm + (times(5) - times(4)) + (times(7) - times(6))
    end do

    ! Dismiss worker processes.

    times(8) = MPI_WTIME()

    time_total = times(8) - times(1)
  else
    ! Worker process.

    times(2) = MPI_WTIME()
    time_setup = times(2) - times(1)

    ! Receive direction if there are/aren't more tasks to perform.
    call MPI_RECV(all_tasks_distributed, 1, MPI_LOGICAL, master_id, tag, &
        MPI_COMM_WORLD, status, err)

    ! Work until no more tasks to be distributed.
    do while (.not. all_tasks_distributed)
      times(3) = MPI_WTIME()

      ! Collect task to complete.
      call MPI_RECV(task, 1, MPI_INTEGER, master_id, tag, MPI_COMM_WORLD, &
          status, err)

      times(4) = MPI_WTIME()

      ! Complete task.
      x_task(:) = 0.0
      call mandelbrot_calculation(N, maxiter, loop_min(task), loop_max(task), &
          x_task)

      times(5) = MPI_WTIME()

      ! Return completed task.
      call MPI_SEND(x_task, chunksize, MPI_REAL, master_id, tag, &
          MPI_COMM_WORLD, err)

      ! Receive direction if there are/aren't more tasks to perform.
      call MPI_RECV(all_tasks_distributed, 1, MPI_LOGICAL, master_id, tag, &
          MPI_COMM_WORLD, status, err)

      times(6) = MPI_WTIME()

      time_wait = time_wait + times(4) - times(3)
      time_comp = time_comp + times(5) - times(4)
      time_comm = time_comm + times(6) - times(5)
    end do

    times(7) = MPI_WTIME()

    time_total = times(7) - times(1)
  end if

  call write_timing_data (N, maxiter, chunksize, n_proc, proc_id, &
      time_setup, time_comp, time_wait, time_comm, time_total)

  ! Writing data to file (only done by master process).
  if (proc_id == master_id) then

    write (*, *) "Writing mandelbrot_master_worker.ppm"

    open(7, file="mandelbrot_master_worker.ppm", status="unknown")

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

  ! MPI finalisation.
  call MPI_FINALIZE(err)

  close(7)

contains

  ! Mandelbrot calculation.
  subroutine mandelbrot_calculation (N, maxiter, lower_bound, upper_bound, x)
    integer , intent(in) :: N, maxiter, lower_bound, upper_bound
    real , intent(out) :: x(0:upper_bound - lower_bound)
    complex :: z, kappa
    integer :: loop, k

    do loop = lower_bound, upper_bound
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

      x(loop-lower_bound) = log(real(k))/log(real(maxiter))
    end do
  end subroutine mandelbrot_calculation

  ! Read in the value of N, maxiter, chunksize from command line arguments.
  subroutine read_input (N, maxiter, chunksize)
    integer , intent(out) :: N, maxiter, chunksize
    integer :: num_args
    character(len=20) :: arg

    num_args = command_argument_count()

    if (num_args < 3) then
      write (*, *) "Usage: <N> <maxiter> <chunksize>"
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

    if (num_args >= 3) then
      call get_command_argument(3, arg)
      read (arg, *) chunksize
    else
      write (*, *) &
          "<chunksize> not specified. Using default value, chunksize = 100000."
      chunksize = 100000
    end if

  end subroutine read_input

  subroutine write_timing_data (N, maxiter, chunksize, n_proc, proc_id, &
      time_setup, time_comp, time_wait, time_comm, time_total)
    integer , intent(in) :: N, maxiter, chunksize, proc_id
    real , intent(in) :: time_setup, time_comp, time_wait, time_comm, time_total
    character(len=1000) :: timing_file
    character(len=20) :: str
    integer :: file_unit

    ! Construct timing filename to be of the form:
    ! "timing-master_worker-N=<N>.maxiter=<maxiter>.n_proc=<n_proc>\
    ! .proc_id=<proc_id>.dat"
    write (timing_file, *) "timing-master_worker-"

    read (str, *) N
    write (timing_file, *) timing_file, "N=", trim(adjustl(str)), "."

    read (str, *) maxiter
    write (timing_file, *) timing_file, "maxiter=", trim(adjustl(str)), "."

    read (str, *) n_proc
    write (timing_file, *) timing_file, "n_proc=", trim(adjustl(str)), "."

    read (str, *) proc_id
    write (timing_file, *) timing_file, "proc_id=", trim(adjustl(str)), ".dat"

    ! Append the chunksize and timing data to the data file
    file_unit = 10 + proc_id

    open (file_unit, file=timing_file, status="append")

    write (file_unit, *) chunksize, time_setup, time_comp, time_wait, time_comm, &
        time_total

    close (file_unit)

  end subroutine write_timing_data

end program mandelbrot_master_worker
