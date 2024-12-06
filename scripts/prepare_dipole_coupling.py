# prepare dipole for QDyn

from os.path import exists

import matplotlib.pyplot as plt
import numpy as np

import qdyn_analyze as qa

########## input ##########
input_energies = ['a.u.', 'eV'][0]
input_coordinates = ['a.u.', 'Angstrom'][0]
input_file = 'PES.dat'
interpolation = ['linear', 'cubic'][0]
convolution = False  # convolution with Gaussian to smooth functions
plot = True
savefig = True
heatmap = True

########## parameters ##########
evtoau = 0.036749405469679
angstoau = 1.8897259886

########## reading input.q ##########
input_file = 'input.q'
if exists(input_file):
    print(qa.read.input_file.__doc__)
    namelist, dynamics, rank, xngrid, yngrid, zngrid, nstates, use_field = qa.read.input_file(input_file)
else:
    print(qa.read.manual_input.__doc__)
    dynamics, rank, xngrid, yngrid, zngrid, nstates, use_field = qa.read.manual_input()

xmin = namelist['general']['xmin']
xmax = namelist['general']['xmax']

########## reading data file ##########
# if rank == 1:
#     if exists(input_file):
#         data = np.genfromtxt(input_file).T
#         if input_energies == 'eV':
#             data[1, :] = data[1, :] * evtoau
#         if input_coordinates == 'Angstrom':
#             data[0, :] = data[0, :] * angstoau
#         data[1, :] = data[1, :] - np.min(data[1, :])
#     else:
#         print('File %s not found. Exiting..' % input_file)
#         exit(1)
#
#     modified = False
#     if data[0, 0] > xmin:
#         print("xmin %.5f is smaller the in the %s (%.5f). Adjust xmin!" % (xmin, input_file, data[0, 0]))
#         xmin = np.round(data[0, 0], 4)
#         modified = True
#
#     if data[0, -1] < xmax:
#         print("xmax %.5f is bigger the in the %s (%.5f). Adjust xmin!" % (xmax, input_file, data[0, -1]))
#         xmax = np.round(data[0, -1], 4)
#         modified = True
#
#     if modified and exists('input.q'):
#         print('Modifying x-range.')
#         namelist['general']['xmin'] = xmin
#         namelist['general']['xmax'] = xmax
#         f90nml.write(namelist, 'newinput.q', force=True)
#         print('Writing new namelist in file newinput.q with corresponding xmin and xmax\n\n', namelist, '\n')
# elif rank == 2:
#     if exists(input_file):
#         data = np.genfromtxt(input_file).T
#         if input_energies == 'eV':
#             data[2, :] = data[2, :] * evtoau
#         if input_coordinates == 'Angstrom':
#             data[0:2, :] = data[0:2, :] * angstoau
#         data[2, :] = data[2, :] - np.min(data[2, :])

########## creating data ##########
# grids
if rank == 1:
    x = np.linspace(xmin, xmax, xngrid)
elif rank == 2:
    X = np.linspace(xmin, xmax, xngrid)
    Y = np.linspace(ymin, ymax, yngrid)
    x = np.zeros(shape=xngrid*yngrid)
    y = np.zeros(shape=xngrid*yngrid)

    k = 0
    for i in range(xngrid):
        for j in range(xngrid):
            # do not change indexing!
            # x goes with the second index (j) because it corresponds to fortran printing
            # otherwise printing will not work
            x[k] = X[j]
            y[k] = Y[i]
            k += 1

# interpolation
# if rank == 1:
#     if interpolation == 'linear':
#         dip = np.interp(x=x, xp=data[0], fp=data[1])
#     elif interpolation == 'cubic':
#         dip = CubicSpline(data[0], data[1], bc_type='natural')(x)
#     else:
#         print("Unknown interpolation technique")
#         exit(1)
# elif rank == 2:
#     if interpolation == 'linear':
#         dip = griddata((data[0], data[1]), data[2], (x, y), method='linear')
#     elif interpolation == 'cubic':
#         dip = griddata((data[0], data[1]), data[2], (x, y), method='cubic')
#     else:
#         print("Unknown interpolation technique")
#         exit(1)

# convolution with Gaussian to smooth out functions
# if convolution:
#     if rank == 1:
#         kernel = np.exp(-(x - np.median(x)) ** 2 / np.var(x) * 1000)
#         kernel = kernel / np.trapz(kernel)
#         dip[10:-10] = np.convolve(dip, kernel, 'same')[10:-10]

# HERE manipulations with potential can be performed
dip = 0.04*np.exp(-x**2/20)
# HERE manipulations with potential can be performed

if plot:
    if rank == 1:
        # plt.plot(data[0], data[1])
        plt.scatter(x, dip, s=5, color='black')
        plt.xlabel('x / a.u.')
        plt.ylabel('E / a.u.')
        if savefig:
            plt.savefig("1D", dpi=300)
        plt.show()
    elif rank == 2:
        if heatmap:
            fig, axs = plt.subplots(1, 1, figsize=(6, 5))
            Z = np.reshape(dip, (xngrid, yngrid))
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
            axs[1].plot_trisurf(x, y, dip, color='C1', alpha=1)
            # axs[1].scatter(x,y,z, color='C1', alpha=1)
            axs[1].set_title('interpolated data')
            plt.tight_layout()
            if savefig:
                plt.savefig("2D", dpi=300)
            plt.show()

# writing pot.dat
if rank == 1:
    np.savetxt('dipole_coup.1.1.dat', np.reshape(dip, newshape=(1, xngrid)), delimiter=' ')
elif rank == 2:
    np.savetxt('dipole_coup.1.1.dat', np.reshape(dip, newshape=(1, xngrid**2)), delimiter=' ')
