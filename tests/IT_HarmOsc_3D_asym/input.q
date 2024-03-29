! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=700,                   ! Number of steps
  dt=1.0,                      ! Timestep [a.u.]
  dtwrite=50.0,                ! Printing every time unit (modulo)

  xngrid=16,                    ! Number of grid points (power of 2 for FFT)
  yngrid=32,                    ! Number of grid points (power of 2 for FFT)
  zngrid=32,                    ! Number of grid points (power of 2 for FFT)
  xmin=-17.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=17.0,
  ymin=-17.0,                  
  ymax=17.0,
  zmin=-17.0,                  
  zmax=17.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_y=2.0,                  ! Reduced mass of system [a.u.]
  mass_z=3.0,                  ! Reduced mass of system [a.u.]
  print_wf=.false.

  rank=3,                      ! Dimensionality
  nstates=10,
/

&it
  pot='0.005*x**2 + 0.0049*y**2 + 0.02535*z**2'  ! Potential
/

&init_wf
/
