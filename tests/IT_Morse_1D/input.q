! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=2500,                  ! Number of steps
  dt=0.5,                      ! Timestep [a.u.]
  dtwrite=10.0,                ! Printing every time unit (modulo)

  xngrid=256,                   ! Number of grid points (power of 2 for FFT)
  xmin=1.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=6.0,
  mass_x=51352.0,                    ! Reduced mass of system [a.u.]
  nstates=10
  print_wf=.false.

  rank=1,                      ! Dimensionality
/

&it
  pot='-0.364 + 0.364 *(1 - exp(-1.51*(-2.041 + x)))^2'  ! Potential
/

&init_wf
/
