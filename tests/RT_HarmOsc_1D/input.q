! Input file for program Qdyn

&general
  run=0,                       ! Type of job (0 - real time propagation, 1 - imaginary time propagation)
  nstep=4000,                  ! Number of steps
  dt=0.1,                      ! Timestep [a.u.]
  dtwrite=1.0,                ! Printing every time unit (modulo)

  ngrid=512,                   ! Number of grid points (power of 2 for FFT)
  xmin=-25.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=25.0,
  mass=5.0,                    ! Reduced mass of system [a.u.]
  wf=0,                        ! Initial wavefunction (0 - generated by program, 1 - wf.chk file)
  nstates=1
  print_wf=.false.
  rank=1,                      ! Dimensionality
  pot='0.005*(x)**2'  ! Potential
/

&init_wf
x0 = 2.0d0
/
