! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=500,                   ! Number of steps
  dt=1.00,                      ! Timestep [a.u.]
  dtwrite=100.0,                 ! Printing every time unit (modulo)

  xngrid=16,                   ! Number of grid points (power of 2 for FFT)
  yngrid=32,                   ! Number of grid points (power of 2 for FFT)
  zngrid=32,                   ! Number of grid points (power of 2 for FFT)
  xmin=-14.0,                  
  xmax=14.0,
  ymin=-13.0,                  
  ymax=13.0,
  zmin=-13.0,                  
  zmax=13.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_y=1.0,                  ! Reduced mass of system [a.u.]
  mass_z=1.0,                  ! Reduced mass of system [a.u.]
  wf=0,                        ! Initial wavefunction (0 - generated by program, 1 - wf.chk file)
  print_wf=.false.

  rank=3,                      ! Dimensionality
  nstates=10,
/

&it
  pot='0.005*x**2 + 0.005*y**2 + 0.005*z**2'  ! Potential
/

&init_wf
/
