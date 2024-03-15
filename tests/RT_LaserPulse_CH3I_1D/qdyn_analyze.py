"""Functions and modules for analysis of QDyn results"""

from os import makedirs, remove, rmdir
from os.path import exists

import f90nml
import imageio
import numpy as np
from matplotlib.pyplot import savefig


# TODO: it would be nice to use as the only input the namelist and the functions would work with just the name list


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
        if rank >= 2:
            namelist['general']['yngrid'] = int(input('yngrid: '))
            namelist['general']['ymin'] = float(input('ymin: '))
            namelist['general']['ymax'] = float(input('ymax: '))
        if rank >= 3:
            namelist['general']['zngrid'] = int(input('zngrid: '))
            namelist['general']['zmin'] = float(input('zmin: '))
            namelist['general']['zmax'] = float(input('zmax: '))
            namelist['rt']['field_coupling'] = bool(input('field_coupling: '))

        return namelist

    def wf(rank, nstates, xngrid, yngrid=0, adiabatic=False):
        wf = []
        for i in range(0, nstates):
            if rank == 1:
                if adiabatic:
                    wf_file = 'wf1d_ad.%d.out'%(i + 1)
                else:
                    wf_file = 'wf1d.%d.out'%(i + 1)
                grid_size = xngrid
            elif rank == 2:
                if adiabatic:
                    wf_file = 'wf2d_ad.%d.out'%(i + 1)
                else:
                    wf_file = 'wf2d.%d.out'%(i + 1)
                grid_size = xngrid*yngrid

            if exists(wf_file):
                wf.append(np.genfromtxt(wf_file))
                # estimate number of frames
                if np.shape(wf[i])[0]%grid_size != 0:
                    print("Number of lines and number of grid point are not matching. Exiting..")
                    exit(1)
                nframes_wf = int(np.shape(wf[i])[0]/grid_size)
                if rank == 1:
                    wf[i] = np.reshape(wf[i], (nframes_wf, grid_size, 5)).transpose((0, 2, 1))
                elif rank == 2:
                    wf[i] = np.reshape(wf[i], (nframes_wf, grid_size, 6)).transpose((0, 2, 1))
            else:
                print("File '%s' does not exist. Exiting.."%wf_file)
                exit(1)

        return np.array(wf), nframes_wf

    def energies(dynamics, nstates):
        en_states = []
        if dynamics == 'rt':
            en_file = 'energies.dat'
            if exists(en_file):
                en_states.append(np.genfromtxt(en_file).T)
                print("'%s' read"%en_file)
            else:
                print("File '%s' does not exist. Exiting.."%en_file)
                exit(1)
            nframes = np.shape(en_states)[2]
        elif dynamics == 'it':
            for i in range(0, nstates):
                en_file = 'energies.%d.dat'%(i + 1)
                if exists(en_file):
                    en_states.append(np.genfromtxt(en_file).T)
                    print("'%s' read"%en_file)
                else:
                    print("File '%s' does not exist. Exiting.."%en_file)
                    exit(1)
            en_states = np.array(en_states)
            print("\nTotal energies:", en_states[:, 1, -1])
            nframes = np.shape(en_states)[2]
        return en_states, nframes

    def field():
        input_file = 'field.dat'
        if not exists(input_file):
            print(f'ERROR: {input_file:s} not found!')
            exit(1)
        return np.genfromtxt(input_file).T

    def pop():
        input_file = 'pop_ad.dat'
        if not exists(input_file):
            print(f'ERROR: {input_file:s} not found!')
            exit(1)
        pop_ad = np.genfromtxt(input_file).T

        input_file = 'pop_diab.dat'
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
        x = np.zeros(shape=xngrid*yngrid)
        y = np.zeros(shape=xngrid*yngrid)

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


class analyze:
    """This class is used to analyze wave functions."""
    def bracket(bra, c, ket):
        # I probably need to convert wf to complex numbers
        return np.trapz(bra*c*ket)


class plot:
    """This class contains functions to plot wave function of its norm."""
