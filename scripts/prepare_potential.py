# prepare potential for QDyn

from os.path import exists

import f90nml
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline, griddata

import qdyn_analyze as qa

########## input ##########
analytic = True  # if True, one needs to set analytic expression below, if Flase, file with pes will be read

if not analytic:
    pes_file = 'PES.dat'
    input_energies = ['a.u.', 'eV'][0]
    input_coordinates = ['a.u.', 'Angstrom'][0]
    interpolation = ['linear', 'cubic'][0]
    shift = True  # shift to minimum of states
    convolution = False  # convolution with Gaussian to smooth functions
    harmonic_approx = True  # creates Harmonic approximation of the potential

plot = True
savefig = False
heatmap = True  # if 2D used

########## parameters ##########
evtoau = 0.036749405469679
angstoau = 1.8897259886

########## reading input.q ##########
input_file = 'input.q'
if exists(input_file):
    print(qa.read.input_file.__doc__)
    namelist = qa.read.input_file(input_file)
else:
    print(qa.read.manual_input.__doc__)
    namelist = qa.read.manual_input()

dynamics = namelist['general']['dynamics']
rank = namelist['general']['rank']
nstates = namelist['general']['nstates']
xngrid = namelist['general']['xngrid']
xmin = namelist['general']['xmin']
xmax = namelist['general']['xmax']
yngrid, zngrid = 0, 0
if rank >= 2:
    yngrid = namelist['general']['yngrid']
    ymin = namelist['general']['ymin']
    ymax = namelist['general']['ymax']
if rank >= 3:
    zngrid = namelist['general']['zngrid']
    zmin = namelist['general']['zmin']
    zmax = namelist['general']['zmax']

########## creating grid based on the input ##########
if rank == 1:
    x = qa.generate.grid_1d(xmin, xmax, xngrid)
elif rank == 2:
    x, y, X, Y = qa.generate.grid_2d(xmin, xmax, xngrid, ymin, ymax, yngrid)

########## Processing input file if not analytic ##########
if not analytic:
    ### reading data file
    if rank == 1:
        if exists(pes_file):
            data = np.genfromtxt(pes_file).T
            if input_energies == 'eV':
                data[1:, :] = data[1:, :]*evtoau
            if input_coordinates == 'Angstrom':
                data[0, :] = data[0, :]*angstoau
            if shift:
                data[1:, :] = data[1:, :] - np.min(data[1, :])
        else:
            print('File %s not found. Exiting..'%pes_file)
            exit(1)

        nstates = np.shape(data)[0] - 1
        print("Number of states: %d"%nstates)
        if nstates > 1:
            state = int(input("State of interest: "))
            if state < 1 or state > nstates:
                print("Only states from 1 to %d are possible!"%nstates)
                exit(1)
        elif nstates == 1:
            state = 1
            print("State of interest: %d"%state)
        else:
            exit(1)

        ### modifying input.q if axes sizes dont match
        modified = False
        if data[0, 0] > xmin:
            print("xmin %.5f is smaller the in the %s (%.5f). Adjust xmin!"%(xmin, pes_file, data[0, 0]))
            xmin = np.round(data[0, 0], 4)
            modified = True

        if data[0, -1] < xmax:
            print("xmax %.5f is bigger the in the %s (%.5f). Adjust xmin!"%(xmax, pes_file, data[0, -1]))
            xmax = np.round(data[0, -1], 4)
            modified = True

        if modified and exists('input.q'):
            print('Modifying x-range.')
            namelist['general']['xmin'] = xmin
            namelist['general']['xmax'] = xmax
            f90nml.write(namelist, 'newinput.q', force=True)
            print('Writing new namelist in file newinput.q with corresponding xmin and xmax\n\n', namelist, '\n')
    elif rank == 2:
        if exists(pes_file):
            data = np.genfromtxt(pes_file).T
            if input_energies == 'eV':
                data[2, :] = data[2, :]*evtoau
            if input_coordinates == 'Angstrom':
                data[0:2, :] = data[0:2, :]*angstoau
            data[2, :] = data[2, :] - np.min(data[2, :])

    ### interpolation
    if rank == 1:
        if interpolation == 'linear':
            pot = np.interp(x=x, xp=data[0], fp=data[state])
        elif interpolation == 'cubic':
            pot = CubicSpline(data[0], data[state], bc_type='natural')(x)
        else:
            print("Unknown interpolation technique")
            exit(1)
    elif rank == 2:
        if interpolation == 'linear':
            pot = griddata((data[0], data[1]), data[2], (x, y), method='linear')
        elif interpolation == 'cubic':
            pot = griddata((data[0], data[1]), data[2], (x, y), method='cubic')
        else:
            print("Unknown interpolation technique")
            exit(1)

    ### convolution with Gaussian to smooth out functions
    if convolution:
        if rank == 1:
            kernel = np.exp(-(x - np.median(x))**2/np.var(x)*1000)
            kernel = kernel/np.trapz(kernel)
            pot[10:-10] = np.convolve(pot, kernel, 'same')[10:-10]
        else:
            print('ERROR: Convolution available only for rank 1!')
            exit(1)

    ### calculating harmonic approximation for the potential
    if harmonic_approx:
        argmin = np.argmin(data[state])
        fit = np.polyfit(data[0, argmin - 1:argmin + 2], data[state, argmin - 1:argmin + 2], 2)
        ho = (fit[0]*x**2 + fit[1]*x + fit[2])
        print(f"Force constant of harmonic oscillator: {fit[0]*2:.6f}")

####### Analytic potential or manipulation with existing potential possible #######
# if not analytic, then you can modify the interpolated potential here
# pot = pot + 0.04*np.exp(20*(x-x[-1]))

### POTENTIAL MODIFICATIONS ###

# if analytic, then this is the part where you should define your potential using coordinate x
# pot = 0.04*np.exp(20*(x-x[-1]))

### ANALYTIC POTENTIAL DEFINITION ###
pot = 0.5*x**2

# note we always work with atomic units


####### Plotting of the potential #######
if plot:
    if rank == 1:
        if not analytic:
            plt.plot(data[0], data[state])
            plt.title("State: %d"%state)
            plt.scatter(x, pot, s=5, color='black')
            fig_name = "state%d_1D"%state

            if harmonic_approx:
                plt.plot(x, ho, linestyle='--', label='HO')
                plt.ylim(np.min(pot), np.max(pot))
        else:
            plt.plot(x, pot, color='black')
            fig_name = 'analytic_pes'

        plt.xlabel('x / a.u.')
        plt.ylabel('E / a.u.')

        if savefig:
            plt.savefig(fig_name, dpi=300)
        plt.show()
    elif rank == 2:
        if heatmap:
            fig, axs = plt.subplots(1, 1, figsize=(6, 5))
            Z = np.reshape(pot, (xngrid, yngrid))
            # pc = axs.pcolormesh(X, Y, np.log(Z+0.05), cmap='hsv')
            pc = axs.pcolormesh(X, Y, Z, cmap='viridis')
            fig.colorbar(pc, ax=axs)
            if savefig:
                plt.savefig("heatmap", dpi=300)
            plt.show()
        else:
            fig = plt.figure(figsize=(10, 6))
            axs = [fig.add_subplot(121, projection='3d'), fig.add_subplot(122, projection='3d')]

            axs[0].plot_trisurf(data[0], data[1], data[2])
            axs[0].set_title('original data')
            axs[1].plot_trisurf(data[0], data[1], data[2], alpha=0.2)
            axs[1].plot_trisurf(x, y, pot, color='C1', alpha=1)
            # axs[1].scatter(x,y,z, color='C1', alpha=1)
            axs[1].set_title('interpolated data')
            plt.tight_layout()
            if savefig:
                plt.savefig("2D", dpi=300)
            plt.show()

# writing pot.dat
if rank == 1:
    if dynamics == 'rt':
        rt_file = 'H.1.1.dat'
        print(f'Saving to file {rt_file:s}')
        np.savetxt(rt_file, np.reshape(pot, newshape=(1, xngrid)), delimiter=' ')
    elif dynamics == 'it':
        np.savetxt('pot.dat', np.reshape(pot, newshape=(1, xngrid)), delimiter=' ')
elif rank == 2:
    if dynamics == 'it':
        np.savetxt('pot.dat', np.reshape(pot, newshape=(1, xngrid*yngrid)), delimiter=' ')
