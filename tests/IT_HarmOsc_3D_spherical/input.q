! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=500,                   ! Number of steps
  dt=1.00,                      ! Timestep [a.u.]
  dtwrite=100.0,                 ! Printing every time unit (modulo)

  xngrid=32,                   ! Number of grid points (power of 2 for FFT)
  yngrid=32,                   ! Number of grid points (power of 2 for FFT)
  zngrid=32,                   ! Number of grid points (power of 2 for FFT)
  xmin=-14.0,                  
  xmax=14.0,
  ymin=-13.0,                  
  ymax=13.0,
  zmin=-13.0,                  
  zmax=13.0,
  mass_x=2.0,                  ! Reduced mass of system [a.u.]
  mass_y=2.0,                  ! Reduced mass of system [a.u.]
  mass_z=2.0,                  ! Reduced mass of system [a.u.]
  print_wf=.false.

  rank=3,                      ! Dimensionality
  nstates=12,
/

&it
  pot='0.01*(x**2 + y**2 + z**2)'  ! Potential
/

&init_wf
/
