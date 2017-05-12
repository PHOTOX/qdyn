! Input file for program Qdyn

&general
  run=0,                      ! Type of job (0 - imaginary time propagation)
  nstep=100,                  ! Number of steps
  dt=1.0,                     ! Timestep [a.u.]

  ngrid=4,                   ! Number of grid points (power of 2 for FFT)
  rank=1,                     ! Dimensionality
  xmin=-1.0,                   ! Grid xmin, xmax same for all dimensions
  xmax=1.0, 
  mass=1.0,                   ! Reduced mass of system [a.u.]

  wf=0,                       ! Initial wavefunction (0 - generated by program, 1 - wf.in file)
  pot='x**2',                ! Potential
/

