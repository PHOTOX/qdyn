# QDYN - quantum dynamics code
Program for one to three-dimensional numerical quantum propagation on a grid in real and imaginary time.

## Compilation
1) Compile FFTW3 library for fast Fourier transform
   - `cd fftw`
   - `tar -xvf fftw.tar.gz` to extract files from the archive, two versions available (3.3.10 and 3.3.6)
   - `./compile.sh` to compile the version selected within the script (3.3.10 default)
   - compilation will output a path of the compilation that will used in Makefile, e.g. `LIBS = -L/home/janos/Programs/fftw/fftw-3.3.6/lib -lfftw3 -lm `
2) Compile Qdym
   - `cd src`
   - change the library path from the previous FFTW3 compilation in `Makefile`
   - `make clean && make` for compilation
   - Ideally, add the `\src` folder to your path so that you can call Qdyn by just `$ qdyn` in the command line
 3) Run tests
    - `make test` - it runs `tests/run_test_suite.sh` which can be accessed separately (python3 with numpy library is necessary for some tests, otherwise they will be skipped)  

To run qdyn, input.q file must be in the folder

## Inputs

Sample input for imaginary time propagation is provided below.
```
&general
  run=1,                       ! Type of job (0 - real time propagation, 1 - imaginary time propagation)
  nstep=200,                   ! Number of steps
  dt=0.5,                      ! Timestep [a.u.]
  dtwrite=10.0,                ! Printing every time unit (modulo)
  ngrid=512,                   ! Number of grid points (power of 2 for FFT)
  xmin=-37.0,                  ! Grid xmin, xmax same for all dimensions
  xmax=37.0,
  mass_x=1.0,                  ! Reduced mass of system [a.u.]
  wf=0,                        ! Initial wavefunction (0 - generated by program, 1 - wf.chk file)
  nstates=10,                  ! Number of states to be optimized
  print_wf=.false.,            ! Printing wavefunction turned off
  rank=1,                      ! Dimensionality
/

&it
  pot='0.005*(x)**2'  ! Potential
  project_rot=.true.
/

&init_wf
/
```

Input is separated into sections:
#### `&general`
The general section sets general variables like time step, number of steps, mass, etc.

#### `&init_wf`
Settings of the initial wave packet.

#### `&it`
Imaginary-time propagation settings.

#### `&rt`
Real-time propagation settings.

## Preparing input
Inputs can be prepared with the help of python scripts.

## Plotting
A series of python scripts is prepared for plotting the data.

## TODO
Imaginary-time propagation and real-time propagation in 1D are finished and tested. Extension to multiple states is currently under progress.
1) Dipole coupling with field
 - read dipoles
2) Transition dipole couplings
3) Diabatic couplings
4) Add autocorrelation function

## small TODO list
1) create ymin, ymax atd.
2) energy functions don't need back FFT
