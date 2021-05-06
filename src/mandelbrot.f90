program mandelbrot

  complex :: z, kappa
  integer :: green, blue, i, j, k, loop

  integer :: N, maxiter
  real , allocatable :: x(:)

  ! Timing variables.
  !!  times       An array storing time markers used to determine the following
  !!              timing variables.
  !!  time_total  Time taken overall.
  double precision :: times(1:2)
  double precision :: time_total

  ! Read command line arguments.
  call read_input(N, maxiter)

  allocate(x(0:N*N-1))

  call cpu_time(times(1))

  ! mandelbrot calculation
  do loop = 0, N*N-1
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

    x(loop) = log(real(k))/log(real(maxiter))
  end do

  call cpu_time(times(2))

  ! Timing analysis.
  time_total = times(2) - times(1)

  write (*, *) &
      "timing for serial code:", NEW_LINE('a'), &
      "  total: ", time_total

  ! writing data to file
  write(*, *) "Writing mandelbrot.ppm"

  open(7, file="mandelbrot.ppm", status="unknown")

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

end program mandelbrot
