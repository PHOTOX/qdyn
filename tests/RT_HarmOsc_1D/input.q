! Input file for program Qdyn

&general
  dynamics='rt',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=4000,                  ! Number of steps
  dt=0.1,                      ! Timestep [a.u.]
  dtwrite=1.0,                ! Printing every time unit (modulo)

  xngrid=512,                   ! Number of grid points (power of 2 for FFT)
  xmin=-10.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=10.0,
  mass_x=10.0,                    ! Reduced mass of system [a.u.]
  nstates=1
  print_wf=.true.
  rank=1,                      ! Dimensionality
/

&init_wf
  x0 = 2.0d0
/

&rt
  pot='0.5*10.0*0.08**2*x**2'  ! Potential 1/2*m*omega^2*x^2
/
