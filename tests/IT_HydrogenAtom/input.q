! Input file for program Qdyn

&general
  dynamics='it',               ! Type of job ('rt' - real time propagation, 'it' - imaginary time propagation)
  nstep=1500,                   ! Number of steps
  dt=0.20,                      ! Timestep [a.u.]
  dtwrite=5.0,                 ! Printing every time unit (modulo)

  xngrid=128,                   ! Number of grid points (power of 2 for FFT)
  yngrid=128,                   ! Number of grid points (power of 2 for FFT)
  zngrid=128,                   ! Number of grid points (power of 2 for FFT)
  xmin=-20.0,                  
  xmax= 20.0,
  ymin=-20.0,                  
  ymax= 20.0,
  zmin=-20.0,                  
  zmax= 20.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  mass_y=1.0,                  ! Reduced mass of system [a.u.]
  mass_z=1.0,                  ! Reduced mass of system [a.u.]
    print_wf=.false.

  rank=3,                      ! Dimensionality
  nstates=5 ,
/

&it
  pot='-1/(x**2 + y**2 + z**2 + 0.000001)^0.5'  ! Potential
/

&init_wf
x0=0.05
y0=0.08
z0=0.10
xsigma=2.0
ysigma=3.0
zsigma=2.5
/
