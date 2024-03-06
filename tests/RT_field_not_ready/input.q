! Input file for program Qdyn

&general
  dynamics='rt',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=4000,                  ! Number of steps
  dt=0.1,                      ! Timestep [a.u.]
  dtwrite=1.0,                ! Printing every time unit (modulo)

  xngrid=512,                   ! Number of grid points (power of 2 for FFT)
  xmin=-25.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=25.0,
  mass_x=5.0,                  ! Reduced mass of system [a.u.]
  nstates=1
  print_wf=.false.
  rank=1,                      ! Dimensionality
/

&init_wf
  x0 = 0.0d0
/

&rt
  pot='0.005*(x)**2'  ! Potential
  field_coupling=.true.
  field='0.001*cos(0.01*t)'
/
