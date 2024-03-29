! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=350,                   ! Number of steps
  dt=0.75,                      ! Timestep [a.u.]
  dtwrite=10.0,                 ! Printing every time unit (modulo)

  xngrid=32,                   ! Number of grid points (power of 2 for FFT)
  yngrid=64,                   ! Number of grid points (power of 2 for FFT)
  xmin=-17.0,                  
  xmax=17.0,
  ymin=-17.0,                  
  ymax=17.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_y=2.0,                  ! Reduced mass of system [a.u.]
  print_wf=.false.

  rank=2,                      ! Dimensionality
  nstates=10
/

&it
  pot='0.005*x**2 + 0.0225*y**2'  ! Potential
/

&init_wf
x0=4.0
y0=4.0
xsigma=0.5
ysigma=0.5
px0=0.0
py0=0.0
/
