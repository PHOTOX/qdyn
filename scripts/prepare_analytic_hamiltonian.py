"""Script for preparing analytic hamiltonian for QDYN calculations.
Qdyn is able to read analytic potentials on for a 1-state system. For multiplte states, input in files H.x.y.dat and
dipole_coup.x.y.dat is required. There is an empty space where the user should define the Hamiltonian and dipole couplings.
Note that the user must ensure that both matrices are symmetric. The script will plot the Hamiltonian, adiabatic energies,
and dipole couplings. The script saves the Hamiltonian and dipole couplings to files readable by QDyn.

© Jiri Janos 2025
"""

import matplotlib.pyplot as plt
import numpy as np

import qdyn_analyze as qa


########## parameters ##########
evtoau = 0.036749405469679
angstoau = 1.8897259886

########## reading input.q ##########
print('\n * Reading input.q')
namelist = qa.read.input_file('input.q')

# creating variables based on the input
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
print('\n * Creating variables')
if rank == 1:
    x = qa.generate.grid_1d(xmin, xmax, xngrid)
elif rank == 2:
    x, y, X, Y = qa.generate.grid_2d(xmin, xmax, xngrid, ymin, ymax, yngrid)
    print("ERROR: rank 2 not implemented yet.")
    exit(1)

# todo: missing rank 2 and 3


### Defining variables ###
if rank == 1:
    H = np.zeros((nstates, nstates, xngrid))  # diabatic Hamiltonian
    E = np.zeros((nstates, xngrid))  # adiabatic energies
    dip = np.zeros((nstates, nstates, xngrid))  # dipole couplings

### ANALYTIC POTENTIAL DEFINITION ###
# here you should define your Hamiltonian and dipole couplings
# e.g. H[0, 0] = 0.10797762*(1 - np.exp(-0.44143377*(x - 0.03984253)))**2
# e.g. H[1, 1] = 0.69575666*np.exp(-0.84386515*(x + 3.05562874)) + 0.12220106  # 3Q0
# e.g. H[2, 2] = 1.2461332*np.exp(-0.59093065*(x + 4.41621395)) + 0.09102091  # 1Q1
# e.g. H[1, 2] = 0.001*np.ones(np.shape(x))
# e.g. dip[0, 1] = 0.1289*np.ones(np.shape(x))
# e.g. dip[0, 2] = 0.1400*np.ones(np.shape(x))
# e.g. dip[1, 2] = 0.1200*np.ones(np.shape(x))

H[0, 0] = 0.10797762*(1 - np.exp(-0.44143377*(x - 0.03984253)))**2
H[1, 1] = 0.69575666*np.exp(-0.84386515*(x + 3.05562874)) + 0.12220106  # 3Q0
H[2, 2] = 1.2461332*np.exp(-0.59093065*(x + 4.41621395)) + 0.09102091  # 1Q1
H[1, 2] = 0.002*np.ones(np.shape(x))
H[2, 1] = H[1, 2]
dip[0, 1] = 0.1289*np.ones(np.shape(x))
dip[1, 0] = dip[0, 1]
dip[0, 2] = 0.1400*np.ones(np.shape(x))
dip[2, 0] = dip[0, 2]
dip[1, 2] = 0.1200*np.ones(np.shape(x))
dip[2, 1] = dip[1, 2]
### ANALYTIC POTENTIAL DEFINITION ###


# checking if matrices are symmetric
print('\n * Checking potentials')
error = False
for i in range(nstates):
    for j in range(i + 1, nstates):
        if any(H[i, j] != H[j, i]):
            print(f'ERROR: H[{i:d}, {j:d}] and H[{j:d}, {i:d}] are not equal! Hamiltonian must be symmetric!')
            error = True
        if any(dip[i, j] != dip[j, i]):
            print(f'ERROR: dip[{i:d}, {j:d}] and dip[{j:d}, {i:d}] are not equal! Dipole couplings must be symmetric!')
            error = True
if error:
    exit(1)

# diagonalization
print('\n * Calculating adiabatic states')
if rank == 1:
    for i in range(xngrid):
        E[:, i] = np.linalg.eigvalsh(H[:, :, i])

### Preparing plot ###
print('\n * Plotting')
plt.rcParams["font.family"] = 'Helvetica'
fig, ax = plt.subplots(3, 1, figsize=(5, 9), sharex='all', gridspec_kw={
    'height_ratios': [1, 1, 0.8]})
state_colors = plt.cm.viridis(np.linspace(0, 1, nstates))
interstate_colors = plt.cm.viridis(np.linspace(0, 1, nstates*(nstates - 1)//2))

### plotting ###
# first plotting diagonal elements
if rank == 1:
    for state in range(nstates):
        ax[0].plot(x, H[state, state], color=state_colors[state], label=f'$H_{{{state:d}{state:d}}}$')
        ax[1].plot(x, E[state], color=state_colors[state], label=f'$E_{{{state:d}}}$')
        ax[2].plot(x, dip[state, state], color=state_colors[state], label=f'$\mu_{{{state:d}{state:d}}}$')
    # then plotting off-diagonal elements
    counter = 0
    for i in range(nstates):
        for j in range(i + 1, nstates):
            ax[0].plot(x, H[i, j], color=interstate_colors[counter], ls='--', label=f'$H_{{{i:d}{j:d}}}$')
            ax[2].plot(x, dip[i, j], color=interstate_colors[counter], ls='--', label=f'$\mu_{{{i:d}{j:d}}}$')
            counter += 1

    ax[0].set_xlim(xmin, xmax)
    ax[0].set_ylabel(r"Diabatic Hamiltonian (a.u.)")
    ax[1].set_xlim(xmin, xmax)
    ax[1].set_ylabel(r"Adiabatic Hamiltonian (a.u.)")
    ax[2].set_xlim(xmin, xmax)
    ax[2].set_xlabel(r"$r$ (a.u.)")
    ax[2].set_ylabel(r"Diab. dipole coupling (a.u.)")

for i in range(3):
    ax[i].minorticks_on()
    ax[i].tick_params('both', direction='in', which='both', top=True, right=True)
    ax[i].legend(frameon=False, labelspacing=0)

plt.tight_layout()
fig.subplots_adjust(hspace=0)
plt.savefig('hamiltonian', dpi=300)
plt.show()

########## saving to files ##########
print('\n * Writing files')
if rank == 1:
    # first save Hamiltonian
    for i in range(nstates):
        for j in range(nstates):
            if all(H[i, j] == 0):  # don't save files if all elements are zero since Qdyn considers that by default
                print(f'H[{i},{j}] = 0, thus, no file is saved.')
                continue
            rt_file = f'H.{i + 1}.{j + 1}.dat'
            print(f'Saving H[{i},{j}] to file {rt_file:s}.')
            np.savetxt(rt_file, np.reshape(H[i, j], shape=(1, xngrid)), delimiter=' ')
    # then save dipole couplings
    for i in range(nstates):
        for j in range(nstates):
            if all(dip[i, j] == 0):  # don't save files if all elements are zero since Qdyn considers that by default
                print(f'dip[{i},{j}] = 0, thus, no file is saved.')
                continue
            rt_file = f'dipole_coup.{i + 1}.{j + 1}.dat'
            print(f'Saving dip[{i},{j}] to file {rt_file:s}.')
            np.savetxt(rt_file, np.reshape(dip[i, j], shape=(1, xngrid)), delimiter=' ')
