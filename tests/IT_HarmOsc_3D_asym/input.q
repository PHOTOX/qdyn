! Input file for program Qdyn

&general
  run=1,                       ! Type of job (0 - real time propagation, 1 - imaginary time propagation)
  nstep=800,                   ! Number of steps
  dt=1.00,                      ! Timestep [a.u.]
  dtwrite=5.0,                 ! Printing every time unit (modulo)

  ngrid=32,                   ! Number of grid points (power of 2 for FFT)
  xmin=-17.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=17.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_y=2.0,                  ! Reduced mass of system [a.u.]
  mass_z=3.0,                  ! Reduced mass of system [a.u.]
  wf=0,                        ! Initial wavefunction (0 - generated by program, 1 - wf.chk file)
  print_wf=.false.

  rank=3,                      ! Dimensionality
  nstates=10,
/

&it
  pot='0.005*x**2 + 0.0049*y**2 + 0.02535*z**2'  ! Potential
/

&init_wf
/
