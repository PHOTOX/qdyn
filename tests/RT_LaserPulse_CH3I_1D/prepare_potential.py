# prepare potential for QDyn

from os.path import exists
import numpy as np

import qdyn_analyze as qa

########## input ##########
analytic = True  # if True, one needs to set analytic expression below, if Flase, file with pes will be read

########## parameters ##########
evtoau = 0.036749405469679
angstoau = 1.8897259886

########## reading input.q ##########
input_file = 'input.q'
if exists(input_file):
    print(qa.read.input_file.__doc__)
    namelist = qa.read.input_file(input_file)
else:
    exit(1)

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
if 'field_coupling' in namelist['rt']:
    use_field = namelist['rt']['field_coupling']
else:
    use_field = False

########## creating grid based on the input ##########
if rank == 1:
    x = qa.generate.grid_1d(xmin, xmax, xngrid)

### ANALYTIC POTENTIAL DEFINITION ###
V_1_1 = 0.10797762*(1 - np.exp(-0.44143377*(x - 0.03984253)))**2
V_2_2 = 0.69575666*np.exp(-0.84386515*(x + 3.05562874)) + 0.12220106
dip = 0.1289*np.ones(np.shape(x))

pot = V_2_2

# writing pot.dat
if rank == 1:
    if dynamics == 'rt':
        rt_file = 'H.1.1.dat'
        print(f'Saving to file {rt_file:s}')
        np.savetxt(rt_file, np.reshape(V_1_1, newshape=(1, xngrid)), delimiter=' ')
        rt_file = 'H.2.2.dat'
        print(f'Saving to file {rt_file:s}')
        np.savetxt(rt_file, np.reshape(V_2_2, newshape=(1, xngrid)), delimiter=' ')
        rt_file = 'dipole_coup.2.1.dat'
        print(f'Saving to file {rt_file:s}')
        np.savetxt(rt_file, np.reshape(dip, newshape=(1, xngrid)), delimiter=' ')
        rt_file = 'dipole_coup.1.2.dat'
        print(f'Saving to file {rt_file:s}')
        np.savetxt(rt_file, np.reshape(dip, newshape=(1, xngrid)), delimiter=' ')
