!> mandelbrot_master_worker
!
!  Adaption of Mandelbrot program to use MPI with master_worker parallelism.
!
!  For a NxN grid (represented in a 1D array; that is, {0, .., N*N - 1}),
!  with n_proc processes, and a specified chunksize.
!
!  The data is broken up into contiguous subsets, with size chunksize, indexed
!  by an integer. Each subset defined by
!  > {loop_min(task), ..., loop_max(task)} for task = 1, ..., n_task.
!
!  The master process will distribute integers to the worker processes,
!  representing which task they are to be working on. The master process will
!  retain a ledger of which task each worker process is working on. If a worker
!  process is not working on any task, this will be recorded in the ledger as
!  task_ledger(proc) = no_task := 0.
!  Initially, all worker processes are assigned a task to work on. The master
!  thread will then wait for the workers to return the Mandelbrot data, and if
!  there are more tasks left to work on, distribute them to the idle workers.
!  When there are no more tasks left to work on, the master process will then
!  switch to collecting unfinished tasks, and telling the workers that return
!  them that there is no more to do. When all outstanding tasks have been
!  collected, the master thread finishes the master-worker scheme.
!
!  The worker processes will first receive a logical flag from the master thread
!  indicating if there is any work to do. If there is, they will then receive an
!  integer indicating which task (that is, which subset of the data) they will
!  work on. After they have performed the Mandelbrot calculation for that subset
!  of data, they will send the Mandelbrot data back to the master thread and
!  wait for its response. When the master thread accepts the workers data, it
!  will tell the worker if there is more work to do, and if there is, send it
!  another task to perform. When there are no more tasks to work on, the worker
!  will finish.
program mandelbrot_master_worker

  use mpi

  complex :: z, kappa
  integer :: green, blue, i, j, k, loop

  integer :: N, maxiter
  real , allocatable :: x(:)

  ! MPI variables.
  !   master_id
  !   proc_id         ID of current MPI process.
  !   n_proc          Number of MPI processes.
  !   err             Stores error code for MPI calls.
  !   tag             MPI tag variable.
  !   request         MPI request variable.
  !   status          MPI status variable. Used by the master process to
  !                   determine the proc_id of worker processes returning tasks.
  !   chunksize       The number of loop iterations assigned to a process for
  !                   one task.
  !   n_tasks         The number of tasks.
  !   loop_min        An array of the lower loop-iteration bounds for each
  !                   task.
  !   loop_max        An array of the upper loop-iteration bounds for each
  !                   task.
  !   x_task          A work array, local to each process, which will yield the
  !                   mandelbrot data for the data subset assigned to this
  !                   worker process. The master thread will use this to copy
  !                   across the mandelbrot data returned by the worker
  !                   processes.
  !   task            A counter variable for looping over tasks.
  !   proc            A counter variable for looping over processes.
  !   proc_recv       The proc_id of the worker process returning a task to the
  !                   master process.
  !   task_recv       The task completed by the worker process returning a task
  !                   to the master process.
  !   task_ledger     A record (kept by the master process) of which task each
  !                   process is currently working on. If the worker process,
  !                   proc, isn't working on any task, then
  !                   task_ledger(proc) = no_task.
  !   no_task         An integer representing an idle task; that is, no task.
  !
  !   all_tasks_distributed   A flag indicating whether or not all the tasks
  !                   have been distributed amongst the worker process. Once the
  !                   master thread has determined that all tasks have been
  !                   distributed, it will send the worker threads this flag as
  !                   they come to return tasks.
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
  !   times       An array storing time markers used to determine the following
  !               timing variables.
  !   time_setup  Time taken for this process to setup MPI variables for
  !               partitioning data.
  !   time_comp   Time taken for a worker process to perform mandelbrot
  !               calculations for its given data subset. For the master thread,
  !               this is used to track the time it takes to copy returned data
  !               into the final data set.
  !   time_wait   For a worker process, this measures the time spent waiting for
  !               the master process to receive this workers completed task.
  !               For the master process, this is not a suitable variable,
  !               since time spent waiting is hard to differentiate from time
  !               spent communicating.
  !   time_comm   For a worker process, this measures the time spent
  !               communicating with the master thread. For the master thread,
  !               this measures the all time spent communicating with the
  !               worker process
  !   time_total  Time taken overall.
  double precision :: times(1:8)
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

        task_ledger(proc) = no_task
      else
        call MPI_ISEND(task, 1, MPI_INTEGER, proc, tag, MPI_COMM_WORLD, &
            request, err)

        task_ledger(proc) = task
        task = task + 1
      end if

      all_tasks_distributed = (task > n_tasks)
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

      times(5) = MPI_WTIME()

      x(loop_min(task_recv):loop_max(task_recv)) = x_task(:)

      times(6) = MPI_WTIME()

      ! Tell worker there is more work to do.
      call MPI_ISEND(all_tasks_distributed, 1, MPI_LOGICAL, proc_recv, tag, &
          MPI_COMM_WORLD, request, err)

      ! Distribute new task.
      call MPI_ISEND(task, 1, MPI_INTEGER, proc_recv, tag, MPI_COMM_WORLD, &
          request, err)

      task_ledger(proc_recv) = task
      task = task + 1

      all_tasks_distributed = (task > n_tasks)

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

      times(5) = MPI_WTIME()

      x(loop_min(task_recv):loop_max(task_recv)) = x_task(:)

      times(6) = MPI_WTIME()

      ! Tell worker there is no more work to do.
      call MPI_ISEND(all_tasks_distributed, 1, MPI_LOGICAL, proc_recv, tag, &
          MPI_COMM_WORLD, request, err)

      ! No more work to distribute.
      task_ledger(proc_recv) = no_task

      times(7) = MPI_WTIME()

      time_comp = time_comp + times(6) - times(5)
      time_comm = time_comm + (times(5) - times(4)) + (times(7) - times(6))
    end do

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

  ! Timing analysis.
  if (proc_id == master_id) then
    write (*, *) &
        "timing for MPI master-worker code, chunksize:", chunksize, &
        NEW_LINE('a'), &
        "  total: ", time_total

    write (*, *) "time spent working/waiting/communicating:"
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, err)

  if (proc_id /= master_id) then
    write (*, '(a, i1, a, f5.2, a, f5.2, a, f5.2, a)') &
        "  ", proc_id, ": ", &
        100.0*time_comp/time_total, " % / ", &
        100.0*time_wait/time_total, " % / ", &
        100.0*time_comm/time_total, " %"
  end if

  ! Write timing data to an output file.
  call write_timing_data (N, maxiter, chunksize, n_proc, proc_id, &
      time_setup, time_comp, time_wait, time_comm, time_total)

  ! Writing Mandelbrot data to file (only done by master process).
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

  ! Write the timing data (for the given parameters: N, maxiter, chunksize,
  ! n_proc), for a given process.
  !
  ! The timing data includes the time taken spent: setting up, communicating,
  ! performing computations, waiting, and the total time spent.
  !
  ! The filename is defined by (N, maxiter, n_proc, proc_id), and the chunksize
  ! is included in the data output to allow for comparison of timing with
  ! varying chunksize.
  subroutine write_timing_data (N, maxiter, chunksize, n_proc, proc_id, &
      time_setup, time_comp, time_wait, time_comm, time_total)
    integer , intent(in) :: N, maxiter, chunksize, proc_id
    double precision , intent(in) :: time_setup, time_comp, time_wait, &
        time_comm, time_total
    character(len=1000) :: timing_file
    character(len=20) :: str_N, str_maxiter, str_n_proc, str_proc_id
    integer :: file_unit

    ! Construct timing filename to be of the form:
    ! "output/timing.master_worker.N=<N>.maxiter=<maxiter>.n_proc=<n_proc>\
    ! .proc_id=<proc_id>.dat"
    write (str_N, *) N
    write (str_maxiter, *) maxiter
    write (str_n_proc, *) n_proc
    write (str_proc_id, *) proc_id

    write (timing_file, *) &
        "output/timing.master_worker.", &
        "N-", trim(adjustl(str_N)), ".", &
        "maxiter-", trim(adjustl(str_maxiter)), ".", &
        "n_proc-", trim(adjustl(str_n_proc)), ".", &
        "proc_id-", trim(adjustl(str_proc_id)), ".dat"

    ! Append the chunksize and timing data to the data file
    file_unit = 10 + proc_id

    open (file_unit, file=trim(adjustl(timing_file)), action="write", &
        position="append")

    write (file_unit, *) chunksize, time_setup, time_comp, time_wait, &
        time_comm, time_total

    close (file_unit)

  end subroutine write_timing_data

end program mandelbrot_master_worker
