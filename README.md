# QDyn: Quantum Dynamics code

**Autors:** Jiří Janoš, Jiří Suchan, Petr Slavíček

Program for one up to three-dimensional numerical quantum propagation on a grid in real and imaginary time.
It is based on the split-operator method, allows for laser pulse interaction and is able to provide exact factorization quantities.

## Features

- Imaginary-time propagation for 1D-3D (allows to optimize arbitrary number of eigenstates).
- Real-time propagation for 1D-3D on a single state (able to calculate the autocorrelation function).
- Real-time propagation in 1D on multiple diabatic states with the laser pulse interaction and exact factorization quantities.

## Compilation

1) Compile FFTW3 library for fast Fourier transform
    - `cd fftw`
    - `tar -xvf fftw.tar.gz` to extract files from the archive, two versions available (3.3.10 and 3.3.6)
    - `./compile.sh` to compile the version selected within the script (3.3.10 default)
    - compilation will output a path of the compilation that will used in Makefile, e.g.
      `LIBS = -L/home/janos/Programs/fftw/fftw-3.3.6/lib -lfftw3 -lm -llapack`
2) Compile QDyn
    - `cd src`
    - change the library path from the previous FFTW3 compilation in `Makefile`
    - `make clean && make` for compilation
    - Ideally, add the `\src` folder to your path so that you can call Qdyn by just `$ qdyn` in the command line
3) Run tests
    - `make test` - it runs `tests/run_test_suite.sh` which can be accessed separately (python3 with numpy library is necessary for
      some tests, otherwise they will be skipped)

To run qdyn, input.q file must be in the folder

## Inputs

Input is provided in `input.q` file, which is read by the program. The input file is divided into sections, each section starts with
`&section_name` and ends with `/`. The sections are described below.
Sample input for real-time dynamics is provided below.

```
&general
    dynamics = 'rt'
    nstep = 350
    dt = 10.0
    dtwrite = 50.0
    xngrid = 2048
    xmin = -18.0
    xmax = 32.0
    mass_x = 2000.0
    nstates = 3
    print_wf = .false.
    rank = 1
/

&init_wf
    x0 = -15.0
    xsigma = 0.7071067811865475
    px0 = 20.0
    init_state = 3
/

&rt
    analytic = .false.
    field_coupling = .false.
/
```

Loads of examples can be found in `tests/` folder.

### `&general`

The `&general` section is mandatory for each calculation and sets general variables like time step, number of steps, mass, etc.

- `dynamics`: Specifies the type of dynamics, e.g. `dynamics='rt'` for real-time dynamics, `dynamics='it'` for imaginary-time
  dynamics.
- `nstep`: Number of time steps for the propagation.
- `dt`: Time step for the propagation in atomic units (a.u.).
- `dtwrite`: Time step for writing the wave function to the output file.
- `rank`: Specifies the rank of the wave function, e.g. `rank=1` for a one-dimensional wave function, `rank=2` for a
  two-dimensional wave function, `rank=3` for a three-dimensional wave function.
- `nstates`: Number of states in the system. For real-time dynamics, this is the number of diabatic states, for imaginary-time
  dynamics, this is the number of eigenstates to be optimized.
- `xngrid`: Number of grid points in the x-direction.
- `yngrid`: Number of grid points in the y-direction (optional, only for 2D and 3D).
- `zngrid`: Number of grid points in the z-direction (optional, only for 3D).
- `xmin`: Minimum value of the x-grid.
- `ymin`: Minimum value of the y-grid (optional, only for 2D and 3D).
- `zmin`: Minimum value of the z-grid (optional, only for 3D).
- `xmax`: Maximum value of the x-grid.
- `ymax`: Maximum value of the y-grid (optional, only for 2D and 3D).
- `zmax`: Maximum value of the z-grid (optional, only for 3D).
- `mass_x`: Mass of the particle in the x-direction in atomic units (a.u.).
- `mass_y`: Mass of the particle in the y-direction in atomic units (a.u.) (optional, only for 2D and 3D).
- `mass_z`: Mass of the particle in the z-direction in atomic units (a.u.) (optional, only for 3D).
- `print_wf`: If set to `.true.`, the wave function will be printed to the output file at each `dtwrite` step.

#### `&init_wf`

Settings of the initial wave packet.

- `gen_init_wf`: If set to `.true.`, the initial wave function will be generated as a Gaussian wave packet. If set to `.false.`,
  the initial wave function will be read from the `wf.chk` file (this file is produced at the end of IT dynamics for the last
  optimized state: it is useful when starting from an eigenstate, e.g., the ground state).
- `init_state`: Specifies in which state should be the wave function guess initialized, e.g. `init_state=2` means that the wave
  function guess will be put in the second state. Applies only to RT dynamics.
- `weights`: Takes the initial wave function guess and spreads across the states with the given weights. Applies only for
  RT dynamics and overrides the use of `init_state`. The input is a list of weights for each state, e.g. `weights=0.1,0.2,0.3` for
  three states. The weights are always renormalized so integer ratios can
  be used also. If the number of weights is less than the number of states, the rest will be set to zero. If the number of weights
  is larger than the number of states, the code will issue error.
- `x0`: Initial position of the Gaussian wave packet in the x-direction.
- `y0`: Initial position of the Gaussian wave packet in the y-direction (optional, only for 2D and 3D).
- `z0`: Initial position of the Gaussian wave packet in the z-direction (optional, only for 3D).
- `xsigma`: Width of the Gaussian wave packet in the x-direction.
- `ysigma`: Width of the Gaussian wave packet in the y-direction (optional, only for 2D and 3D).
- `zsigma`: Width of the Gaussian wave packet in the z-direction (optional, only for 3D).
- `px0`: Initial momentum of the Gaussian wave packet in the x-direction.
- `py0`: Initial momentum of the Gaussian wave packet in the y-direction (optional, only for 2D and 3D).
- `pz0`: Initial momentum of the Gaussian wave packet in the z-direction (optional, only for 3D).

#### `&it`

Imaginary-time propagation settings.

- `analytic`: If set to `.true.`, the code will use an analytic potential specified in the `pot` keyword. If set to `.false.`,
  code will read Hamiltonian from `pot.dat` file (can be produced with `scripts/prepare_analytic_hamiltonian.py`).
- `pot`: Specifies the potential to be used in the imaginary-time propagation if `analytic` is set to `.true.`. It should be a
  string representing the potential, e.g. `pot='x^2'` for a harmonic oscillator potential.
- `project_rot`: If set to `.true.`, the code projects out also lower eigenstate rotated by 90 degrees in the complex plane. If
  set to `.false.`, the code will only project out only the lower eigenstates which may lead to optimizing the same eigenstate
  twice just rotate by 90 degrees. (default `.true.`)

#### `&rt`

Real-time propagation settings.

- `analytic`: If set to `.true.`, the code will use an analytic potential specified in the `pot` keyword (only for 1D dynamics). If
  set to `.false.`, code will read Hamiltonian from `H.x.x.dat` file (can be produced with
  `scripts/prepare_analytic_hamiltonian.py`).
- `pot`: Specifies the potential to be used in the real-time propagation if `analytic` is set to `.true.`. It should be a
  string representing the potential, e.g. `pot='x^2'` for a harmonic oscillator potential.
- `field_coupling`: If set to `.true.`, the code will use the field coupling term in the Hamiltonian.
- `field`: Specifies the field to be used in the real-time propagation if `field_coupling` is set to `.true.`. It should be a
  string representing the field, e.g. `field='sin(0.1*t)'` for a sinusoidal field.
- `field_on`: Specifies the time when the field is turned on, e.g. `field_on=0.0` means the field is turned on at the beginning of
  the propagation.
- `field_off`: Specifies the time when the field is turned off, e.g. `field_off=100.0` means the field is turned off after 100 a.u.
  of time.
- `norm_thresh`: Threshold for the normalization of the wave function, e.g. `norm_thresh=1e-6` means that QDyn stops when the
  norm starts differing from 1 by more than `norm_thresh`.
- `exact_factor`: If set to `.true.`, the code will calculate exact factorization quantities. (available only for 1D RT dynamics)
- `ef_gauge`: Specifies the gauge for the exact factorization quantities, `ef_gauge='A0'` for the $\vec{A}=0$ gauge,
  `ef_gauge='S0'` for the $S=0$ gauge.
- `ef_zero`: Effective EF zero, i.e., when density is below this value, the code will not calculate the exact factorization
  quantities.
  This is useful for avoiding numerical noise in the calculation of the exact factorization quantities.
- `autocorr`: If set to `.true.`, the code will calculate the autocorrelation function, outputed to `autocorr.dat` file.

## Preparing QDyn inputs and processing of data

QDyn is supplemented with a set of Python scripts that can be used to prepare the input files and process the output data.
The key script is `qdyn_analyze.py` which contains all functions for preparing inputs and processing the output data from QDyn.
There are functions to read all the output properties and calculate several quantities. The remaining scripts are specific
scripts drawing funcitons from `qdyn_analyze.py` to plot the dynamics and calculate data.

## References

This software uses the FFTW3 library (Frigo and Johnson, Proc. IEEE, 93(2), 2005). See http://www.fftw.org.

## Citation
If you use QDyn in your research, please cite:
```bibtex
@software{qdyn,
  title = {QDyn: Quantum Dynamics code},
  author = {Janoš, Jiří; Suchan, Jiří; Slavíček, Petr},
  url = {https://github.com/PHOTOX/qdyn},
  year = {2025},
}
```