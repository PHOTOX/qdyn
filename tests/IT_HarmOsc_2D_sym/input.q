! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=500,                   ! Number of steps
  dt=1.0,                      ! Timestep [a.u.]
  dtwrite=10.0,                 ! Printing every time unit (modulo)

  xngrid=64,                   ! Number of grid points (power of 2 for FFT)
  yngrid=32,                   ! Number of grid points (power of 2 for FFT)
  xmin=-17.0,                  
  xmax=17.0,
  ymin=-14.0,                  
  ymax=16.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_y=1.0,                  ! Reduced mass of system [a.u.]
  print_wf=.false.

  rank=2,                      ! Dimensionality
  nstates=10
/

&it
  pot='0.005*x**2 + 0.005*y**2'  ! Potential
/

&init_wf
/
