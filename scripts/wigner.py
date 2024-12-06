"""Perform wigner transformation of the wave function"""

from os.path import exists

import matplotlib.pyplot as plt
import numpy as np

import qdyn_analyze as qa

########## INPUT ##########
wig_state = 0  # state for wigner transformation, 0 = ground state
dy = 0.005  # numerical integration step for evaluating the Wigner transformation integral
plot_density = True  # plot Wigner density function
nsamples = 0  # number of samples taken from the Wigner density, if nsamples<=0, then no sampling is performed

# defaults
input_file = 'input.q'

########## reading input.q ##########
if exists(input_file):
    print(qa.read.input_file.__doc__)
    namelist = qa.read.input_file(input_file)
else:
    print(qa.read.manual_input.__doc__)
    namelist = qa.read.manual_input()

# saving input data to variables
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

########## initialize ##########
if dynamics == 'rt':
    print("ERROR: Real time dynamics not allowed for wigner tranformation.")
    exit(1)

if rank != 1:
    print("ERROR: Currently implemented only for 1 dimension")
    exit(1)

# getting data
wf, pot_en, [x, y, z], nframes_wf = qa.read.wf(rank=rank, nstates=nstates, xngrid=xngrid, yngrid=yngrid)

wigner_wf = wf[wig_state, -1]  # wave function to be wigner transformed (last step of optimization for state=wig_state)

if plot_density:
    # calculating wigner transformed density wigner transformation
    wigner_transform = qa.calc.wigner1d_density(wigner_wf, x, dy=dy)

    # plotting wigner transformed density
    plt.rcParams["font.family"] = 'Helvetica'
    fig, axs = plt.subplots(1, 1, figsize=(6, 5))

    xw, pw = np.meshgrid(np.unique(wigner_transform[0]), np.unique(wigner_transform[1]))
    zw = np.reshape(wigner_transform[2], np.shape(xw), order='F')

    pc = axs.pcolormesh(xw, pw, zw, cmap='RdBu')
    cm_max = np.max(np.abs(zw))  # setting max and min for colorbar so that 0 is white, red negative and blue positive
    pc.set_clim(-cm_max, cm_max)
    fig.colorbar(pc, ax=axs)

    axs.set_ylabel(r'p (a.u.)')
    axs.set_xlabel(r'x (a.u.)')
    axs.set_title('Wigner density')
    plt.tight_layout()
    plt.show()

if nsamples > 0:
    # IMPORTANT: this sampling generates x and p from the grid plotted above. This won't be effective as the grid is
    # very wide so the user is advised to use optional keywords [xmin, xmax, pmin, pmax] to make the range narrower and
    # sampling more effective. These can be selected based on the plot of distributions.

    # selecting range from which the random numbers [0, rnd_max] for comparing with the density will be drawn
    if 'wigner_transform' in locals():  # if I have the wigner density calculated, I will use the maximum of the desity
        rnd_max = 1.1*np.max(wigner_transform[2])
    else:
        rnd_max = 1/np.pi  # this should be the upper boundary for any Wigner function

    samples = qa.calc.sample_from_wig(wigner_wf, x, nsamples=nsamples, rnd_max=rnd_max, dy=dy)

    # plotting the sampled density and true density in positions
    plt.rcParams["font.family"] = 'Helvetica'
    fig, axs = plt.subplots(1, 1, figsize=(4, 4))

    axs.plot(x, np.conjugate(wigner_wf)*wigner_wf, label=r'$|\psi|^2$')
    axs.hist(samples[0], density=True, alpha=0.2, align='mid', bins=max([int(nsamples/100), 10]), label="sampling")

    axs.set_ylabel(r'density')
    axs.set_xlabel(r'x (a.u.)')
    axs.legend(labelspacing=0.1, frameon=False)
    plt.tight_layout()
    plt.show()

    print("\nWigner sampled positions: \n", samples[0])
    print("\nWigner sampled momenta: \n", samples[1])