# analyze wf from QDyn

import f90nml
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists
from scipy.interpolate import CubicSpline

########## input ##########
input_energies = ['a.u.', 'eV'][0]
input_coordinates = ['a.u.', 'Angstrom'][1]
input_file = 'PES.dat'
interpolation = ['linear', 'cubic spline'][1]
convolution = False # convolution with Gaussian to smooth functions
plot = True

########## parameters ##########
evtoau = 0.036749405469679
angstoau = 1.8897259886

########## reading input.q ##########
if exists('input.q'):
    namelist = f90nml.read('input.q')
    print('\nReading namelist in file input.q\n\n', namelist, '\n')
    ngrid = namelist['general']['ngrid']
    xmin = namelist['general']['xmin']
    xmax = namelist['general']['xmax']
else:
    print('File input.q not found. Using defaults. Modify variables manually in the code.')
    ngrid = int(input('ngrid: '))
    xmin = float(input('xmin: '))
    xmax = float(input('xmax: '))

########## reading data file ##########
if exists(input_file):
    data = np.genfromtxt(input_file).T
    if input_energies == 'eV':
        data[1, :] = data[1, :] * evtoau
    if input_coordinates == 'Angstrom':
        data[0, :] = data[0, :] * angstoau
    data[1, :] = data[1, :] - np.min(data[1, :])
else:
    print('File %s not found. Exiting..' % input_file)
    exit(1)

modified = False
if data[0, 0] > xmin:
    print("xmin %.5f is smaller the in the %s (%.5f). Adjust xmin!" % (xmin, input_file, data[0, 0]))
    xmin = np.round(data[0, 0], 4)
    modified = True

if data[0, -1] < xmax:
    print("xmax %.5f is bigger the in the %s (%.5f). Adjust xmin!" % (xmax, input_file, data[0, -1]))
    xmax = np.round(data[0, -1], 4)
    modified = True

if modified and exists('input.q'):
    print('Modifying x-range.')
    namelist['general']['xmin'] = xmin
    namelist['general']['xmax'] = xmax
    f90nml.write(namelist, 'newinput.q', force=True)
    print('Writing new namelist in file newinput.q with corresponding xmin and xmax\n\n', namelist, '\n')

########## creating data ##########
# xrange created
x = np.linspace(xmin, xmax, ngrid)

# interpolation
if interpolation == 'linear':
    pot = np.interp(x=x, xp=data[0], fp=data[1])
elif interpolation == 'cubic spline':
    pot =  CubicSpline(data[0], data[1], bc_type='natural')(x)
else:
    print("Unknown interpolation technique")
    exit(1)

# convolution with Gaussian to smooth out functions
if convolution:
    kernel=np.exp(-(x - np.median(x)) ** 2 / np.var(x) * 1000)
    kernel= kernel / np.trapz(kernel)
    pot[10:-10]=np.convolve(pot, kernel, 'same')[10:-10]

# HERE manipulations with potential can be performed
pot = pot + 0.04*np.exp(20*(x-x[-1]))
# HERE manipulations with potential can be performed

if plot:
    plt.plot(data[0], data[1])
    plt.scatter(x, pot, s=5, color='black')
    plt.xlabel('x / a.u.')
    plt.ylabel('E / a.u.')
    plt.show()

# noinspection PyTypeChecker
np.savetxt('pot.dat', np.reshape(pot, newshape=(1, ngrid)), delimiter=' ')
