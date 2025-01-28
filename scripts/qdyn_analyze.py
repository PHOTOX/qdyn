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
                    wf_file = f'{folder:s}/wf1d_ad.{i + 1:d}.out'
                else:
                    wf_file = f'{folder:s}/wf1d.{i + 1:d}.out'
                grid_size = xngrid
            elif rank == 2:
                if adiabatic:
                    wf_file = f'{folder:s}/wf2d_ad.{i + 1:d}.out'
                else:
                    wf_file = f'{folder:s}/wf2d.{i + 1:d}.out'
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

    def make_gif(gif_frame, duration=2, gif_folder='gif_frames'):
        imageio.mimsave('wavepacket.gif', gif_frame, duration=duration)
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
        dy = step for numerical integration."""

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

    def quantum_momenta_1d(wf, x, positions):
        "Calculating quantum momenta (Curchod, Martinez, 2018, eq. 37)."

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
            quantum_momenta[i] = np.imag(np.conjugate(wf_interp(xp)) * wf_num_deriv(xp)) / (
                    np.conjugate(wf_interp(xp)) * wf_interp(xp))

        return np.array([positions, quantum_momenta])

    def bracket_1d(bra, c, ket, x):
        """1D bracket calculation."""
        bracket = np.conjugate(bra) * c * ket
        return np.trapz(bracket, x=x)

    def tbf_1d(x, x0, p, alpha, gamma=1 + 0j, amplitude=1 + 0j):
        """Creating a TBF as used in AIMS. It contains amplitude, phase, momentum and prefactors."""
        phase = np.exp(gamma * 1j)
        prefactor = (2 * alpha / np.pi) ** (1 / 4)
        position = np.exp(-alpha * (x - x0) ** 2)
        momentum = np.exp(p * 1j * (x - x0))

        tbf_wf = amplitude * phase * prefactor * position * momentum
        return tbf_wf


def progress(percent, width, n, str=''):
    """Function to print progress of calculation."""
    left = width * percent // n
    right = width - left
    print(f'\r{str}[', '#' * left, ' ' * right, '] %d' % (percent * 100 / n) + '%', sep='', end='', flush=True)
