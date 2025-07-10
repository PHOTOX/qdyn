"""Functions and modules for analysis of QDyn results.

© Jiri Janos
"""

from os import makedirs, remove, rmdir
from os.path import exists

import f90nml
import imageio
import numpy as np
from matplotlib.pyplot import savefig


# TODO: it would be nice to use as the only input the namelist and the functions would work with just the name list
class units:
    """Conversion factors.
    Useful source: https://onlinelibrary.wiley.com/doi/pdf/10.1002/3527605606.app9"""
    # fs to atomic time units
    fstoau = 41.341374575751
    autofs = 1 / fstoau
    # eV to atomic units (Hartree)
    evtoau = 0.036749405469679
    autoev = 1 / evtoau
    # Angstrom to atomic units (Bohr)
    angstoau = 1.8897259886
    autoangs = 1 / angstoau
    # cm^-1 to atomic energy units (Hartree)
    cm1toau = 0.0000045563352812122295
    autocm1 = 1 / cm1toau
    # nanometer to atomic energy units (Hartree)
    nmtoau = 45.56335
    autonm = 1 / nmtoau
    # V/cm to atomic units
    autovcm = 5.14220826e9
    vcmtoau = 1 / autovcm
    # atomic mass unit to atomic units
    amutoau = 1822.8884902498949
    autoamu = 1 / amutoau


class read:
    """Class containing all the functions used for reading data."""

    def input_file(input_file):
        """Reading input file with f90nml package to set qdyn variables.
        dynamics = 'rt' or 'it' # type of dynamics
        rank = 1, 2, 3 # dimension
        xngrid = grid size
        yngrid = grid size
        zngrid = grid size
        nastates >= 1
        use_field = True or False
        """
        if not exists(input_file):
            print(f'ERROR: {input_file:s} not found!')
            exit(1)
        namelist = f90nml.read(input_file)
        print('\nReading namelist in file input.q\n\n', namelist)

        return namelist

    def manual_input():
        """Asking for manual input of qdyn variables.
        dynamics = 'rt' or 'it' # type of dynamics
        rank = 1, 2, 3 # dimension
        xngrid = grid size
        yngrid = grid size
        zngrid = grid size
        nastates >= 1
        use_field = True or False
        """

        # probably not the most fancy way of creating dictionary but there is no connection on the plane
        namelist = dict()
        namelist['general'] = dict()
        namelist['rt'] = dict()
        namelist['it'] = dict()

        while True:
            dynamics = str(input('dynamics: '))
            if dynamics == 'it' or dynamics == 'rt':
                namelist['general']['dynamics'] = dynamics
                break
            else:
                print("Allowed options for dynamics are 'it' and 'rt'.")

        namelist['general']['rank'] = int(input('rank: '))
        namelist['general']['nstates'] = int(input('nstates: '))

        namelist['general']['xngrid'] = int(input('xngrid: '))
        namelist['general']['xmin'] = float(input('xmin: '))
        namelist['general']['xmax'] = float(input('xmax: '))
        if namelist['general']['rank'] >= 2:
            namelist['general']['yngrid'] = int(input('yngrid: '))
            namelist['general']['ymin'] = float(input('ymin: '))
            namelist['general']['ymax'] = float(input('ymax: '))
        if namelist['general']['rank'] >= 3:
            namelist['general']['zngrid'] = int(input('zngrid: '))
            namelist['general']['zmin'] = float(input('zmin: '))
            namelist['general']['zmax'] = float(input('zmax: '))
            namelist['rt']['field_coupling'] = bool(input('field_coupling: '))

        return namelist

    def wf(rank, nstates, xngrid, yngrid=0, adiabatic=False, folder='.'):
        """Reading wave function from the wf files.
        Currently wf_data is old outputing of the wf file as it is read while wf is the complex array."""
        wf_data = []
        for i in range(0, nstates):
            if rank == 1:
                if adiabatic:
                    wf_file = f'{folder:s}/wf1d_ad.{i + 1:d}.dat'
                else:
                    wf_file = f'{folder:s}/wf1d.{i + 1:d}.dat'
                grid_size = xngrid
            elif rank == 2:
                if adiabatic:
                    wf_file = f'{folder:s}/wf2d_ad.{i + 1:d}.dat'
                else:
                    wf_file = f'{folder:s}/wf2d.{i + 1:d}.dat'
                grid_size = xngrid * yngrid

            if exists(wf_file):
                wf_data.append(np.genfromtxt(wf_file))
                # estimate number of frames
                if np.shape(wf_data[i])[0] % grid_size != 0:
                    print("Number of lines and number of grid point are not matching. Exiting..")
                    exit(1)
                nframes_wf = int(np.shape(wf_data[i])[0] / grid_size)
                if rank == 1:
                    wf_data[i] = np.reshape(wf_data[i], (nframes_wf, grid_size, 5)).transpose((0, 2, 1))
                elif rank == 2:
                    wf_data[i] = np.reshape(wf_data[i], (nframes_wf, grid_size, 6)).transpose((0, 2, 1))
                print(f"'{wf_file:s}' read")
            else:
                print("File '%s' does not exist. Exiting.." % wf_file)
                exit(1)

        wf_data = np.array(wf_data)

        # here extract essential data and create a complex wave function
        pot_en = np.zeros((nstates, grid_size))
        wf = np.zeros((nstates, nframes_wf, grid_size), dtype=complex)
        x, y, z = np.zeros((grid_size)), np.zeros((grid_size)), np.zeros((grid_size))

        x = wf_data[0, 0, 0]

        for state in range(nstates):
            pot_en[state] = wf_data[state, 0, -1]
            for frame in range(nframes_wf):
                wf[state, frame] = wf_data[state, frame, -4] + 1j * wf_data[state, frame, -3]

        if rank == 2:
            y = wf_data[0, 0, 1]
        elif rank == 2:
            z = wf_data[0, 0, 2]

        grid = [x, y, z]

        return wf, pot_en, grid, nframes_wf

    def energies(dynamics, nstates):
        en_states = []
        if dynamics == 'rt':
            en_file = 'energies.dat'
            if exists(en_file):
                en_states.append(np.genfromtxt(en_file).T)
                print("'%s' read" % en_file)
            else:
                print("File '%s' does not exist. Exiting.." % en_file)
                exit(1)
        elif dynamics == 'it':
            for i in range(0, nstates):
                en_file = 'energies.%d.dat' % (i + 1)
                if exists(en_file):
                    en_states.append(np.genfromtxt(en_file).T)
                    print("'%s' read" % en_file)
                else:
                    print("File '%s' does not exist. Exiting.." % en_file)
                    exit(1)
        en_states = np.array(en_states)
        if dynamics == 'it': print("\nTotal energies:", en_states[:, 1, -1])
        nframes = np.shape(en_states)[2]
        return en_states, nframes

    def field(input_file='field.dat'):
        if not exists(input_file):
            print(f'ERROR: {input_file:s} not found!')
            exit(1)
        return np.genfromtxt(input_file).T

    def pop(folder='./'):
        input_file = f'{folder:s}pop_ad.dat'
        if not exists(input_file):
            print(f'ERROR: {input_file:s} not found!')
            exit(1)
        pop_ad = np.genfromtxt(input_file).T

        input_file = f'{folder:s}pop_diab.dat'
        if not exists(input_file):
            print(f'ERROR: {input_file:s} not found!')
            exit(1)
        pop_diab = np.genfromtxt(input_file).T

        return pop_ad, pop_diab

    def dipole_coupling(rank, nstates, xngrid, yngrid=0, folder='.'):
        """Reading dipole moments from the dipole_coup.n.m.dat files."""

        if rank == 1:
            dipoles = np.zeros(shape=(nstates, nstates, xngrid))
            for i in range(0, nstates):
                for j in range(0, nstates):
                    dip_file = f'{folder:s}/dipole_coup.{i + 1:d}.{j + 1:d}.dat'
                    if exists(dip_file):
                        dipoles[i, j] = np.genfromtxt(dip_file)
                        print(f"'{dip_file:s}' read")
                    else:
                        print(f"'{dip_file:s}' does not exist. Considering it all zeros.")
        else:
            print("ERROR: Reading dipole moments is currently implemented only for 1 dimension.")
            exit(1)

        return dipoles

    def ef(rank, nstates, xngrid, yngrid=0, folder='.'):
        """Reading Exact Factorization quantities."""

        # get grid size
        if rank == 1:
            grid_size = xngrid
        elif rank == 2:
            grid_size = xngrid * yngrid

        # nuclear density
        nucdens_file = f'{folder:s}/nuclear_density.dat'
        if exists(nucdens_file):
            nucdens = np.genfromtxt(nucdens_file)

            # estimate number of frames
            if np.shape(nucdens)[0] % grid_size != 0:
                print(f"Number of lines and number of grid point are not matching in {nucdens_file:s}. Exiting..")
                exit(1)
            nframes = int(np.shape(nucdens)[0] / grid_size)

            if rank == 1:
                nucdens = np.reshape(nucdens, (nframes, grid_size, 2)).transpose((0, 2, 1))
                x, y = nucdens[0, 0], np.zeros(shape=(xngrid))
                nucdens = nucdens[:, 1]
            elif rank == 2:
                nucdens = np.reshape(nucdens, (nframes, grid_size, 3)).transpose((0, 2, 1))
                x, y = nucdens[0, 0], nucdens[0, 1]
                nucdens = nucdens[:, 2]
            print(f"'{nucdens_file:s}' read")
        else:
            print("File '%s' does not exist. Exiting.." % nucdens_file)
            exit(1)

        grid = [x, y]

        # nuclear phase and gradient
        nucphase_file = f'{folder:s}/nuclear_phase.dat'
        if exists(nucphase_file):
            nucphase = np.genfromtxt(nucphase_file)

            if np.shape(nucphase)[0] % grid_size != 0:
                print(f"Number of lines and number of grid point are not matching in {nucphase_file:s}. Exiting..")
                exit(1)
            elif nframes != int(np.shape(nucphase)[0] / grid_size):
                print(f"Number of time frames in {nucdens_file:s} is not consistent with {nucphase_file:s}. Exiting..")
                exit(1)

            if rank == 1:
                nucphase = np.reshape(nucphase, (nframes, grid_size, 3)).transpose((0, 2, 1))
                nucphase = nucphase[:, 1:]
            elif rank == 2:
                nucphase = np.reshape(nucphase, (nframes, grid_size, 4)).transpose((0, 2, 1))
                nucphase = nucphase[:, 2:]
            print(f"'{nucphase_file:s}' read")
        else:
            print("File '%s' does not exist. Exiting.." % nucphase_file)
            exit(1)

        # GI-TDPES
        gitdpes_file = f'{folder:s}/gi-tdpes.dat'
        if exists(gitdpes_file):
            gitdpes = np.genfromtxt(gitdpes_file)

            if np.shape(gitdpes)[0] % grid_size != 0:
                print(f"Number of lines and number of grid point are not matching in {gitdpes_file:s}. Exiting..")
                exit(1)
            elif nframes != int(np.shape(gitdpes)[0] / grid_size):
                print(f"Number of time frames in {nucdens_file:s} is not consistent with {gitdpes_file:s}. Exiting..")
                exit(1)

            if rank == 1:
                gitdpes = np.reshape(gitdpes, (nframes, grid_size, 6)).transpose((0, 2, 1))
                gitdpes = gitdpes[:, 1:]
            elif rank == 2:
                gitdpes = np.reshape(gitdpes, (nframes, grid_size, 7)).transpose((0, 2, 1))
                gitdpes = gitdpes[:, 2:]
            print(f"'{gitdpes_file:s}' read")
        else:
            print("File '%s' does not exist. Exiting.." % gitdpes_file)
            exit(1)

        # TDVP
        tdvp_file = f'{folder:s}/tdvp.dat'
        if exists(tdvp_file):
            tdvp = np.genfromtxt(tdvp_file)

            # estimate number of frames
            if np.shape(tdvp)[0] % grid_size != 0:
                print(f"Number of lines and number of grid point are not matching in {tdvp_file:s}. Exiting..")
                exit(1)
            nframes = int(np.shape(tdvp)[0] / grid_size)

            if rank == 1:
                tdvp = np.reshape(tdvp, (nframes, grid_size, 2)).transpose((0, 2, 1))
                tdvp = tdvp[:, 1]
            elif rank == 2:
                tdvp = np.reshape(tdvp, (nframes, grid_size, 3)).transpose((0, 2, 1))
                tdvp = tdvp[:, 2]
            print(f"'{tdvp_file:s}' read")
        else:
            print("File '%s' does not exist. Exiting.." % tdvp_file)
            exit(1)

        # electronic coefficients
        C_file = f'{folder:s}/el_coefficients_ef.dat'
        if exists(C_file):
            C = np.genfromtxt(C_file)

            # estimate number of frames
            if np.shape(C)[0] % grid_size != 0:
                print(f"Number of lines and number of grid point are not matching in {C_file:s}. Exiting..")
                exit(1)
            nframes = int(np.shape(C)[0] / grid_size)

            if rank == 1:
                C = np.reshape(C, (nframes, grid_size, 1 + nstates * 2)).transpose((0, 2, 1))
                C = C[:, 1:]
            elif rank == 2:
                C = np.reshape(C, (nframes, grid_size, + nstates * 2)).transpose((0, 2, 1))
                C = C[:, 2:]
            el_coeff = np.zeros(shape=(nframes, nstates, grid_size), dtype=complex)
            for step in range(nframes):
                for state in range(nstates):
                    el_coeff[step, state] = C[step, 2 * state] + 1j * C[step, 2 * state + 1]
            print(f"'{C_file:s}' read")
        else:
            print("File '%s' does not exist. Exiting.." % C_file)
            exit(1)

        # GD-TDPES
        gdtdpes_file = f'{folder:s}/gd-tdpes.dat'
        if exists(gdtdpes_file):
            gdtdpes = np.genfromtxt(gdtdpes_file)

            if np.shape(gdtdpes)[0] % grid_size != 0:
                print(f"Number of lines and number of grid point are not matching in {gdtdpes_file:s}. Exiting..")
                exit(1)
            elif nframes != int(np.shape(gdtdpes)[0] / grid_size):
                print(f"Number of time frames in {nucdens_file:s} is not consistent with {gdtdpes_file:s}. Exiting..")
                exit(1)

            if rank == 1:
                gdtdpes = np.reshape(gdtdpes, (nframes, grid_size, 2)).transpose((0, 2, 1))
                gdtdpes = gdtdpes[:, 1]
            elif rank == 2:
                gdtdpes = np.reshape(gdtdpes, (nframes, grid_size, 3)).transpose((0, 2, 1))
                gdtdpes = gdtdpes[:, 2]
            print(f"'{gdtdpes_file:s}' read")
        else:
            print("File '%s' does not exist. Exiting.." % gdtdpes_file)
            exit(1)

        return nucdens, nucphase, gitdpes, gdtdpes, tdvp, el_coeff, grid


class gif:
    """Class of functions for creating gif. 
    It uses following variables:
    - gif_frames: list storing the figures for gif
    - gif_folder: stores the figs when saved (default available)
    - duration: what's the time gap between frames in gif (default available)
    - gif_dpi: dpi for gif resolution (default available)
    """

    def init_gif(gif_folder='gif_frames'):
        if exists(gif_folder):
            rmdir(gif_folder)
        makedirs(gif_folder)
        gif_frames = []
        return gif_frames

    def save_frame(gif_frames, dpi_gif=100, gif_folder='gif_frames'):
        fig_file = gif_folder + f'/frame.png'
        savefig(fig_file, format='png', dpi=dpi_gif)
        gif_frames.append(imageio.v2.imread(fig_file))
        remove(fig_file)

    def make_gif(gif_frame, duration=2, gif_folder='gif_frames', gif_name='wavepacket'):
        imageio.mimsave(gif_name + '.gif', gif_frame, duration=duration)
        rmdir(gif_folder)


class generate:
    """This class is used to generate data based on input, e.g. grids"""

    def grid_1d(xmin, xmax, xngrid):
        x = np.linspace(xmin, xmax, xngrid)
        return x

    def grid_2d(xmin, xmax, xngrid, ymin, ymax, yngrid):
        X = np.linspace(xmin, xmax, xngrid)
        Y = np.linspace(ymin, ymax, yngrid)
        x = np.zeros(shape=xngrid * yngrid)
        y = np.zeros(shape=xngrid * yngrid)

        k = 0
        for i in range(xngrid):
            for j in range(yngrid):
                # do not change indexing!
                # x goes with the second index (j) because it corresponds to fortran printing
                # otherwise printing will not work
                x[k] = X[j]
                y[k] = Y[i]
                k += 1
        return x, y, X, Y


class calc:
    """This class is used to calculations like expectation values or Wigner sampling"""

    def bracket_1d(bra, c, ket, x):
        """1D bracket calculation."""
        bracket = np.conjugate(bra) * c * ket
        return np.trapz(bracket, x=x)

    def wf_energy_1d(wf, pot, mass, x):
        """Wave function energy. This function is calculates energy of the wave function on a given state in BH
        expansion, it does not calculate energy of the total wave function."""

        # check a single wf, not an array from BH expansion
        if np.ndim(wf) != 1:
            print("ERROR: dimension is not one, probably more than one wf insterted!")
            exit(1)

        # check that lenght of wf is the same as of x
        if np.shape(wf)[0] != np.shape(x)[0]:
            print("ERROR: Length of wf and x are not the same!")
            exit(1)

        if np.shape(wf)[0] != np.shape(pot)[0]:
            print("ERROR: Length of wf and pot are not the same!")
            exit(1)

        # norm
        norm = np.real(calc.bracket_1d(wf, 1, wf, x))

        # potential energy
        V = np.real(calc.bracket_1d(wf, pot, wf, x)) / norm

        # kinetic energy
        p = 2 * np.pi * np.fft.fftfreq(np.shape(x)[0], (x[1] - x[0]))  # momentum space
        wfp = np.fft.ifft(wf)  # Fourier transform of the wave function
        wfx = np.fft.fft(wfp * p ** 2 / 2 / mass)  # kinetic energy operator in momentum space and transform to coordinate space
        T = np.real(calc.bracket_1d(wf, 1, wfx, x)) / norm  # kinetic energy

        # total energy
        energy = V + T

        return energy

    def wigner1d_integral(wf, x, xp, pp, dy):
        """Wigner integral for 1D wave function at point [xp, pp].
        dy is the step for numerical integration."""

        def wf_interp(arg):
            """Interpolated wave function for wigner evaluation"""
            return np.interp(x=arg, xp=x, fp=wf)

        # getting information about the x array
        xmin, xmax, xngrid = np.min(x), np.max(x), np.shape(x)[0]
        # setting y for integration and its limits
        # the limits are bigger than the x size to ensure the integral is almost from -infinity to infinity
        ymin = xmin - 1 * (xmax - xmin)
        ymax = xmax + 1 * (xmax - xmin)
        y = np.arange(ymin, ymax, dy)

        integrand = wf_interp(xp + y / 2) * np.conjugate(wf_interp(xp - y / 2)) * np.exp(pp * y * 1j)
        return np.real(np.trapz(x=y, y=integrand) / 2 / np.pi)

    def wigner1d_density(wf, x, dy=0.002):
        """Calculating wigner transform of a given wave function.
        The calculation considers atomic units. Therefore, h=1 in the Wigner transform.
        dy = step for numerical integration.
        It calls wigner1d_integral function to calculate the integral at each point of the grid."""

        # preparing p range for the Wigner density using the FFT grid
        # the FFT grid should be sufficient as it's used to propagate wf with the split operator technique
        xngrid = np.shape(x)[0]
        p = np.sort(np.fft.fftfreq(xngrid, (x[1] - x[0])) * 2)  # momentum space

        # wigner distribution
        wigner = []
        for xp in x:
            for pp in p:
                wig = calc.wigner1d_integral(wf, x, xp, pp, dy)
                wigner.append([xp, pp, wig])
            progress(np.where(xp == x)[0][0], 50, len(x), str='Generating Wigner distribution: ')
        print("\nDone!")  # printing an empty line after progress

        wigner = np.array(wigner).T

        return wigner

    def sample_from_wig(wf, x, nsamples, dy=0.002, xmin=0, xmax=0, pmin=0, pmax=0, rnd_max=1.0):
        """Sampling from the Wigner distribution.
        dy = step for numeric integration of the Wigner integral
        nsamples = number of samples selected from the distribution
        xmin, xmax, pmin, pmax = max and min values for which we sample; if set to 0, range determined automatically
         from grids that are used during the pslit operator propagation (for p we use the same FFT grid)."""

        # setting xmin, xmax, pmin, pmax if not specified by user (same as when generating the density)
        if xmin == 0:
            xmin = np.min(x)
        if xmax == 0:
            xmax = np.max(x)
        if pmin == 0:  # minimum of FFT grid used in the split operator for propagation
            pmin = np.min(np.fft.fftfreq(np.shape(x)[0], x[1] - x[0]))
        if pmax == 0:  # maximum of FFT grid used in the split operator for propagation
            pmax = np.max(np.fft.fftfreq(np.shape(x)[0], x[1] - x[0]))

        # checking rnd_max
        if rnd_max <= 0:
            print("ERROR: rnd_max is <=0 and cannot be used to select random number from range [0, rnd_max]."
                  "\nPlease use rnd_max keyword to change it.")

        # generating lists of sampled data
        samples = np.zeros(shape=(2, nsamples))
        # drawing the samples
        nattempts = 0
        for i in range(nsamples):
            while True:
                nattempts += 1
                xp = np.random.uniform(low=xmin, high=xmax, size=1)[0]  # random x
                pp = np.random.uniform(low=pmin, high=pmax, size=1)[0]  # random p
                rnd = np.random.uniform(low=0, high=rnd_max, size=1)[0]  # random number to be compared with Wig. dist.

                wig = calc.wigner1d_integral(wf, x, xp, pp, dy)

                # check if the distribution is not negative or higher than rnd_max
                if wig < -1e-4 * rnd_max:  # this threshold is for numerical accuracy
                    print(f"ERROR: Sampling from negative Wigner distribution is not possible!\n"
                          f"For x={xp} and p={pp} the Wigner probability is {wig}.")
                    exit(1)
                elif wig > rnd_max:
                    print(f"ERROR: rnd_max ({rnd_max}) is smaller than Wigner probability ({wig} for x={xp} and "
                          f"p={pp}). Increase rnd_max and rerun again.")
                    exit(1)

                # check if the point is sampled
                if rnd <= wig:
                    samples[0, i] = xp
                    samples[1, i] = pp
                    break

            progress(i + 1, 50, nsamples, str='Sampling from Wigner distribution: ')
        print("\nDone!")  # printing an empty line after progress
        print(f'Success rate of random sampling: {nsamples / nattempts * 100:.6f}%')

        return samples

    def sample_from_density(density, x, nsamples=100, rank=1, xmin=0, xmax=0):
        """Sampling from an inputted density.
        nsamples = number of samples selected from the distribution
        xmin, xmax = max and min values for which we sample; if set to 0, range determined automatically
         from the grid."""

        # check dimension
        if rank != 1:
            print("ERROR: Sampling from density is implemented only for 1 dimension.")
            exit(1)

        # setting xmin, xmax if not specified by user
        if xmin == 0:
            xmin = np.min(x)
        if xmax == 0:
            xmax = np.max(x)

        # setting rnd_max
        rnd_max = np.max(density)

        # generating lists of sampled data
        samples = np.zeros(shape=nsamples)
        # drawing the samples
        for i in range(nsamples):
            while True:
                xp = np.random.uniform(low=xmin, high=xmax, size=1)
                den_xp = np.interp(xp, x, density)
                rnd = np.random.uniform(low=0, high=rnd_max, size=1)
                if rnd < den_xp:
                    samples[i] = xp
                    break

        print(f' * Density sampled!')
        return samples

    def quantum_momenta_1d(wf, x, positions):
        """Calculating quantum momenta (Curchod, Martinez, 2018, eq. 37).
        Alternative definition of the quantum momenta is in Quantum Chemistry and Dynamics of Excited States, 2021,
        Chapter 18, Bohmian Approaches to Non-adiabatic Molecular Dynamics, eq. 18.5.
        Currently, the definition from the book (eq. 18.5) is used as it is thriftier."""

        def wf_interp(arg):
            """Interpolated wave function for wigner evaluation"""
            return np.interp(x=arg, xp=x, fp=wf)

        def wf_num_deriv(arg, dx=0.02):
            """Numerical first derivative."""
            return (wf_interp(arg + dx) - wf_interp(arg)) / dx

        nsamples = np.shape(positions)[0]
        quantum_momenta = np.zeros(shape=np.shape(positions))
        for i in range(nsamples):
            xp = positions[i]
            # Eq 18.5 from Quantum Chemistry and Dynamics of Excited States
            quantum_momenta[i] = np.imag(wf_num_deriv(xp) / wf_interp(xp))
            # alternatively Eq 37 from Curchod, Martinez, 2018, CURRENTLY PHASED OUT SINCE THE OTHER FORMULA IS THRIFTIER
            # quantum_momenta[i] = np.imag(np.conjugate(wf_interp(xp))*wf_num_deriv(xp))/(
            #             np.conjugate(wf_interp(xp))*wf_interp(xp))

        return np.array([positions, quantum_momenta])

    def tbf_1d(x, x0, p, alpha, gamma=1 + 0j, amplitude=1 + 0j):
        """Creating a TBF as used in AIMS. It contains amplitude, phase, momentum and prefactors."""
        phase = np.exp(gamma * 1j)
        prefactor = (2 * alpha / np.pi) ** (1 / 4)
        position = np.exp(-alpha * (x - x0) ** 2)
        momentum = np.exp(p * 1j * (x - x0))

        tbf_wf = amplitude * phase * prefactor * position * momentum
        return tbf_wf

    def el_density_matrix_1d(wf, x):
        """Calculating electron density matrix from the wave function.
        The electronic density is defined as rho_ij = <chi_i|chi_j>.
        Diagonals should be the electronic state populations while off-diagonals are coherences."""

        nstates, ntimes, nxgrid = np.shape(wf)

        if nxgrid != np.shape(x)[0]:
            print("ERROR: The grid size and the wave function size do not match!")
            exit(1)

        density_matrix = np.zeros(shape=(nstates, nstates, ntimes), dtype=complex)

        for t in range(ntimes):
            for i in range(nstates):
                for j in range(nstates):
                    density_matrix[i, j, t] = np.trapz(np.conjugate(wf[i, t]) * wf[j, t], x)

        return density_matrix

    def dipoles_1d(wf, dip_coup, x):
        """Calculating time-dependent expectation value of the dipole moment.
        The expectation value of dipole moment is defined as mu = sum_ij <chi_i|mu_ij|chi_j>.
        :param wf: wave function
        :param dip_coup: dipole coupling matrix
        :param x: grid
        """

        nstates, ntimes, nxgrid = np.shape(wf)

        if nxgrid != np.shape(x)[0] or nxgrid != np.shape(dip_coup)[2]:
            print("ERROR: The grid size and the wave function size do not match!")
            exit(1)

        dipole = np.zeros(shape=(ntimes), dtype=float)

        for t in range(ntimes):
            for i in range(nstates):
                for j in range(i, nstates):
                    if i == j:
                        dipole[t] += np.real(calc.bracket_1d(wf[i, t], dip_coup[i, j], wf[j, t], x))
                    else:
                        dipole[t] += 2 * np.real(calc.bracket_1d(wf[i, t], dip_coup[i, j], wf[j, t], x))

        return dipole

    def mean_x_1d(wf, x):
        """Calculating the expectation value of the position."""
        nstates, ntimes, nxgrid = np.shape(wf)
        mean_x = np.zeros(shape=(ntimes))

        for t in range(ntimes):
            for i in range(nstates):
                mean_x[t] += calc.bracket_1d(wf[i, t], x, wf[i, t], x).real

        return mean_x

    def aver_force_1D(wf, pot_en, dip_coup, x, field, time_unit='a.u.'):
        """Calculating the expectation value of force in 1D as defined in Cardosa-Guitierrez, Remacle, J Phys B 2024.
        The force is decomposed into four terms, F_V, F_[tau,V], F_mu, F_[tau,mu]. Currently, only the terms
        without commutators with NAC (tau) are implemented.
        :param wf: wave function
        :param pot_en: potential energy
        :param dip_coup: dipole coupling matrix
        :param x: grid
        """
        import matplotlib.pyplot as plt

        if time_unit not in ['a.u.', 'fs']:
            print("ERROR: Time unit not recognized. Exiting..")
            exit(1)

        nstates, ntimes, nxgrid = np.shape(wf)
        force_pot, force_dip = np.zeros(shape=(nstates, ntimes), dtype=complex), np.zeros(shape=(ntimes), dtype=complex)

        for t in range(ntimes):
            # force coming from potential
            for i in range(nstates):
                force_pot[i, t] += - calc.bracket_1d(wf[i, t], np.gradient(pot_en[i], x, edge_order=2), wf[i, t], x)

            for i in range(nstates):
                for j in range(nstates):
                    force_dip[t] += calc.bracket_1d(wf[i, t], np.gradient(dip_coup[i, j], x, edge_order=2), wf[j, t], x) * \
                                    field[1, t]

        if np.max(force_dip.imag) > 1e-15 or np.max(force_pot.imag) > 1e-15:
            print("ERROR: Imaginary part of the force is not negligible. Exitting..")
            exit(1)
        else:
            force_dip = np.array(force_dip.real, dtype=float)
            force_pot = np.array(force_pot.real, dtype=float)

        # force_total = force_pot + force_dip
        force_total = np.sum(force_pot, axis=0) + force_dip

        # Calculating the total force as a time derivative of the wave function from the definition (this is what
        # I should get from force_total). This is simpler and better way of getting the total force, however,
        # it doesn't provide the decomposition into different terms.
        if time_unit == 'fs':
            time = field[0] * units.fstoau
        else:
            time = field[0]

        p = np.zeros(shape=(ntimes))
        for t in range(ntimes):
            for i in range(nstates):
                p[t] += calc.bracket_1d(wf[i, t], 1, -1j * np.gradient(wf[i, t], x, edge_order=2), x).real
        force_p = np.gradient(p, time, edge_order=2)

        # check if my total force is correct
        if np.max(np.abs(force_total - force_p)) > 1e-10:
            print(
                "WARNING: The total force is not the same as the time derivative of the wave function.")  # exit(1) # todo: turn it on later

        # todo: the following section is for testing purposes and should be removed later
        # problem at the moment is that my total force calculated from F_V and F_mu is not the same as the time derivative,
        # i.e., the exact formula. Testing showed that the F_V should be correct, there is just small difference which was
        # attributed to the numeric derivative accuracy (I copied the exact derivative formula and it worked better).
        # However, the F_mu term show too large oscillations. I tested by comparing <x> and values obtained by propagating Newton equations.
        # When I use force_p, it all works perfectly. Thus, it seems I have a problem with the F_mu term, which I can't solve right now.
        # Or I possibly need other terms but they should be all zero since they depend on the NAC, which is 0.
        # Thus, I do a dirty trick and get force_mu = force_p - force_pot. This is not ideal but it should work for now.
        testing = False
        if testing:
            # propagate x=1,8 with newton equations using force_total as force
            mean_x = calc.mean_x_1d(wf, x)
            x = np.zeros(shape=(ntimes))
            v = np.zeros(shape=(ntimes))
            x[0] = mean_x[0]

            dt = time[1] - time[0]
            mass = 1728.256714
            for i in range(1, ntimes):
                x[i] = x[i - 1] + v[i - 1] * dt
                v[i] = v[i - 1] + force_total[
                    i - 1] / mass * dt  # v[i] = v[i - 1] + force_p[i - 1]/mass*dt  # v[i] = v[i - 1] + np.sum(force_pot, axis=0)[i - 1]/mass*dt

            plt.subplot(121)
            # plt.plot(time,force_pot[0])
            # plt.plot(time,force_pot[1])
            # plt.plot(time, np.sum(force_pot, axis=0))
            # plt.plot(time, force_dip/field[1])
            plt.plot(time, force_dip, label='calculated $F_{\mu}$')
            # plt.plot(time, force_total)
            # plt.plot(time, np.convolve(force_total, np.ones(20)/20, mode='same'))
            # plt.plot(time, np.convolve(force_p, np.ones(20)/20, mode='same'), ls='--')
            # plt.plot(time, force_p, ls='--')
            plt.plot(time, force_p - np.sum(force_pot, axis=0), ls='-.', label='"correct" $F_{\mu}$')
            # plt.plot(time, force_p - np.sum(force_pot, axis=0) - force_dip, ls='-.')
            # plt.plot(time, (force_p - np.sum(force_pot, axis=0) - force_dip)/(force_dip/field[1]), ls='--')

            plt.legend(labelspacing=0, frameon=False)

            plt.subplot(122)

            plt.plot(time, mean_x)
            # plt.plot(time,p)
            plt.plot(time, x)
            plt.legend(labelspacing=0, frameon=False)
            plt.show()
            plt.pause(100)
            exit()

        print("WARNING: This is only a temporary code and the solution is not guaranteed to be correct.")
        force_total = force_p
        force_dip = force_total - np.sum(force_pot, axis=0)
        # I could possibly make one big array with all the forces and return it
        return force_total, force_pot, force_dip


def progress(percent, width, n, str=''):
    """Function to print progress of calculation."""
    left = width * percent // n
    right = width - left
    print(f'\r{str}[', '#' * left, ' ' * right, '] %d' % (percent * 100 / n) + '%', sep='', end='', flush=True)
