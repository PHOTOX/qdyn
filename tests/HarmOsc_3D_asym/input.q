! Input file for program Qdyn

&general
  run=1,                       ! Type of job (0 - real time propagation, 1 - imaginary time propagation)
  nstep=1500,                   ! Number of steps
  dt=0.50,                      ! Timestep [a.u.]
  dtwrite=5.0,                 ! Printing every time unit (modulo)

  ngrid=64,                   ! Number of grid points (power of 2 for FFT)
  xmin=-15.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=15.0,
  mass=1.0,                    ! Reduced mass of system [a.u.]
  wf=0,                        ! Initial wavefunction (0 - generated by program, 1 - wf.chk file)
  print_wf=.false.

  rank=3,                      ! Dimensionality
  nstates=10,
  pot='0.005*x**2 + 0.00245*y**2 + 0.00845*z**2'  ! Potential

/

&init_wf
/
