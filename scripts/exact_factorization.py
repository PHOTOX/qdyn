"""Extracting exact factorization quantities from QDyn simulations.
The code is currently suitable only for 1D (RT-)dynamics with various gauges.

© Jiri Janos
"""

from os.path import exists

import matplotlib.pyplot as plt
import numpy as np

import qdyn_analyze as qa

########## INPUT ##########
# general plotting settings
pause = 0.00001  # pause between frames during plotting
frame_step = 1  # plot only frame_step instead of plotting every frame
wf_scaling = 0.12  # scaling factor for the wf and density so that they are visible in the plots
plot_online = True  # if true, plots the results on the fly, otherwise it's just the final plot
adiabatic = True  # adiabatic vs diabatic BH states
plot_density = True  # if true, plot density instead of wave function, false plot wave function
xplotrange = [1, 3.5]  # if empty array, it takes xmin and xmax from the input file

# gif
gif = True  # if true, saves the frames for gif
dpi_gif = 100
duration = 2  # in ms, (duration=1000 * 1/fps, fps=frames per second).

# units for plotting
time_unit = ['a.u.', 'fs'][1]  # time units in femtoseconds or atomic time units


########## functions ##########
def plot_ef_1d():
    """Plotting the exact factorization results in 1D."""
    # setting canvas
    fig, axs = plt.subplots(3, 3, figsize=(13, 9), gridspec_kw={'height_ratios': [4, 4, 1]})
    axs_bh = axs[0, 0]
    axs_ef = axs[1, 0]
    axs_field = axs[2, 0]
    axs_S = axs[0, 1]
    axs_gdtdpes = axs[1, 1]
    axs_pop = axs[2, 1]
    axs_elcoef = axs[0, 2]
    axs_tdpes = axs[1, 2]
    axs_en = axs[2, 2]

    # define colors
    colors = plt.cm.viridis([0.0, 0.4, 0.8, 0.2, 0.6, 1.0])

    # getting axis ranges
    if len(xplotrange) == 2:
        xminplot = xplotrange[0]
        xmaxplot = xplotrange[1]
    else:
        xminplot = xmin
        xmaxplot = xmax
    if plot_density:
        vmax = 3 * np.max(den_bh) + np.min(pot_en[-1])
        vmin = np.min([np.min(den_bh), np.min(pot_en)])
    else:
        vmax = 1.5 * np.max(wf_bh)
        vmin = np.min([np.min(wf_bh), np.min(pot_en)])

    for i in range(0, nframes, frame_step):
        # cleaning canvas
        axs_bh.cla()
        axs_ef.cla()
        axs_field.cla()
        axs_S.cla()
        axs_tdpes.cla()
        axs_gdtdpes.cla()
        axs_elcoef.cla()
        axs_pop.cla()
        axs_en.cla()

        fig.suptitle(f'Time: {t[i]:.0f} {time_unit:s}  Energy: {energy[0, i]:.4f} a.u.')

        # plotting Born-Huang quantities
        for j in range(0, nstates):
            axs_bh.plot(x, pot_en[j], color='black')

            if plot_density:
                axs_bh.plot(x, den_bh[j][i] + pot_en[j], color=colors[j], label=rf'$|\psi_{j:d}|^2$')
                axs_bh.fill_between(x, den_bh[j][i] * 0 + pot_en[j], den_bh[j][i] + pot_en[j], color=colors[j],
                                    alpha=0.35)
            else:
                axs_bh.plot(x, np.real(wf_bh[j][i]) + pot_en[j], label=rf'Re[$\psi_{j:d}$]', color=colors[j])
                axs_bh.plot(x, np.imag(wf_bh[j][i]) + pot_en[j], label=rf'Im[$\psi_{j:d}$]', color=colors[j], ls='--')

        axs_bh.set_xlim(xminplot, xmaxplot)
        axs_bh.set_ylim(vmin, vmax)
        axs_bh.set_ylabel(r'$E$ (a.u.)')
        axs_bh.set_xlabel(r'$x$ (a.u.)')
        axs_bh.legend(labelspacing=0)

        # plotting EF quantities
        for j in range(0, nstates):  # BH states
            axs_ef.plot(x, pot_en[j], color='black', alpha=0.2)

        if plot_density:
            # axs_ef.plot(x, nucdens[i] + gitdpes[i, -1], color='C2', label=r'$|\chi|^2 + \varepsilon_{GI}$')
            # axs_ef.fill_between(x, nucdens[i] + gitdpes[i, -1], gitdpes[i, -1], color='C2', alpha=0.35)
            axs_ef.plot(x, nucdens[i] + energy[0, i], color=colors[0], label=r'$|\chi|^2$')
            axs_ef.fill_between(x, nucdens[i] + energy[0, i], energy[0, i], color=colors[0], alpha=0.35)
        else:
            nucwf = np.sqrt(nucdens[i]) * np.exp(1j * nucphase[i, 0])
            axs_ef.plot(x, np.real(nucwf) + energy[0, i], color=colors[0], label=r'Re$[\chi]$')
            axs_ef.plot(x, np.imag(nucwf) + energy[0, i], color=colors[0], ls='--', label=r'Im$[\chi]$')

        axs_ef.plot(x, gitdpes[i, -1], color=colors[1], label=r'$\varepsilon_{GI}$')

        axs_ef.set_xlim(xminplot, xmaxplot)
        axs_ef.set_ylim(vmin, vmax)
        axs_ef.set_ylabel(r'$E$ (a.u.)')
        axs_ef.set_xlabel(r'$x$ (a.u.)')
        axs_ef.legend(labelspacing=0)

        # plotting nuclear phase of TDVP based on the gauge
        if ef_gauge == 'A0':
            # nuclear phase
            axs_S.plot(x, nucphase[i, 0], linewidth=2, color=colors[0], label=f'$S$')
            axs_S.plot(x, nucphase[i, 1], linewidth=1, color=colors[1], label=r'$\nabla S$')
            axs_S.axhline(0, linewidth=0.5, color='black')
            axs_S.set_xlim(xminplot, xmaxplot)
            axs_S.set_ylim(np.min(nucphase[:, 0]), np.max(nucphase[:, 0]))
            axs_S.set_ylabel(r'$S$ (a.u.)')
            axs_S.set_xlabel(r'$x$ (a.u.)')
            axs_S.legend(labelspacing=0)
        elif ef_gauge == 'S0':
            # TDVP
            axs_S.plot(x, tdvp[i], color=colors[0], label=r'$\vec{A}$')
            axs_S.axhline(0, linewidth=0.5, color='black')
            axs_S.set_xlim(xminplot, xmaxplot)
            axs_S.set_ylim(np.min(tdvp[:i + 1]), np.max(tdvp[:i + 1]))
            axs_S.set_ylabel(r'$\vec{A}$ (a.u.)')
            axs_S.set_xlabel(r'$x$ (a.u.)')
            axs_S.legend(labelspacing=0)

        # plotting GI-TDPES
        axs_tdpes.plot(x, gitdpes[i, 4], lw=2, ls='-', color=colors[0], label=r'$\varepsilon_{GI}$', zorder=-1)
        axs_tdpes.plot(x, gitdpes[i, 0], lw=1, ls='-', color=colors[1], label=r'$H_{el}$')
        axs_tdpes.plot(x, gitdpes[i, 1], lw=1, ls='-', color=colors[2], label=r'$V_{int}$')
        axs_tdpes.plot(x, gitdpes[i, 2], lw=1, ls='--', color=colors[5], label=r'$|\nabla C|^2$')
        axs_tdpes.plot(x, gitdpes[i, 3], lw=1, ls='--', color=colors[4], label=r'$A^2$')
        axs_tdpes.axhline(0, linewidth=0.5, color='black')
        axs_tdpes.set_xlim(xminplot, xmaxplot)
        axs_tdpes.set_ylim(np.min(gitdpes[:i + 1]))
        axs_tdpes.set_ylabel(r'$\varepsilon_{GI}$ (a.u.)')
        axs_tdpes.set_xlabel(r'$x$ (a.u.)')
        axs_tdpes.legend(labelspacing=0)

        # GD-TDPES
        axs_gdtdpes.plot(x, gdtdpes[i], color=colors[0], lw=2, label=r'$\varepsilon_{GD}$', zorder=0)
        tdpes = gdtdpes[i] + gitdpes[i, 4]
        axs_gdtdpes.plot(x, tdpes, color=colors[1], lw=2, label=r'$\varepsilon$', zorder=-1)
        averageE = np.trapz(nucdens[i] * tdpes, x)
        if plot_density:
            axs_gdtdpes.plot(x, nucdens[i] + averageE, color=colors[2],
                             label=r'$|\chi|^2 + \langle\chi|\varepsilon|\chi\rangle$')
            axs_gdtdpes.fill_between(x, nucdens[i] + averageE, averageE, color=colors[2], alpha=0.35)
        else:
            nucwf = np.sqrt(nucdens[i]) * np.exp(1j * nucphase[i, 0])
            axs_gdtdpes.plot(x, np.real(nucwf) + averageE, color=colors[2],
                             label=r'Re$[\chi] + \langle\chi|\varepsilon|\chi\rangle$')
            axs_gdtdpes.plot(x, np.imag(nucwf) + averageE, color=colors[2], ls='--',
                             label=r'Im$[\chi] + \langle\chi|\varepsilon|\chi\rangle$')
        axs_gdtdpes.axhline(0, linewidth=0.5, color='black')
        axs_gdtdpes.set_xlim(xminplot, xmaxplot)
        axs_gdtdpes.set_ylim(min([np.min(tdpes), np.min(gdtdpes)]) - 0.2 * np.max(nucdens),
                             max([np.max(nucdens + averageE), np.max(gdtdpes)]) + 0.2 * np.max(nucdens))
        axs_gdtdpes.set_ylabel(r'$\varepsilon_{GD}$ (a.u.)')
        axs_gdtdpes.set_xlabel(r'$x$ (a.u.)')
        axs_gdtdpes.legend(labelspacing=0)

        # el. coefficients
        for j in range(0, nstates):
            axs_elcoef.plot(x, np.abs(el_coeff[i, j]) ** 2, color=colors[j], label=rf'$|C_{j:d}|^2$')
        axs_elcoef.axhline(0, lw=0.5, color='black')
        axs_elcoef.set_xlim(xminplot, xmaxplot)
        axs_elcoef.set_ylim(-0.05, 1.05)
        axs_elcoef.set_ylabel(r'Electronic coefficients (a.u.)')
        axs_elcoef.set_xlabel(r'$x$ (a.u.)')
        axs_elcoef.legend(labelspacing=0)

        # populations
        for state in range(1, nstates + 1):
            axs_pop.plot(t[:i + 1], pop_ad[state][:i + 1], color=colors[state - 1], label=rf'$P_{state - 1:d}$')
        axs_pop.set_ylim(-0.1, 1.1)
        axs_pop.set_xlim(t[0], t[-1])
        axs_pop.set_xlabel(f'$t$ ({time_unit:s})')
        axs_pop.set_ylabel(f'Adiab. pop.')

        # energies
        axs_en.plot(t[:i + 1], energy[0][:i + 1], color=colors[0], label=r'$E_\mathregular{tot}$')
        axs_en.set_xlim(t[0], t[-1])
        xde = (np.max(energy[0]) - np.min(energy[0]))
        axs_en.set_ylim(np.min(energy[0]) - 0.1 * xde, np.max(energy[0]) + 0.1 * xde)
        axs_en.set_xlabel(f'$t$ ({time_unit:s})')
        axs_en.set_ylabel(f'$E$ (a.u.)')
        axs_en.legend()

        # field
        if use_field:
            axs_field.plot(field[0, :i], field[1, :i], color=colors[0])
            axs_field.set_xlim(field[0, 0], field[0, -1])
            axs_field.set_xlim(t[0], t[-1])
            axs_field.set_xlabel(f'$t$ ({time_unit:s})')
            axs_field.set_ylabel(f'el. field (a.u.)')
        else:
            axs_field.plot(t[:i + 1], energy[0][:i + 1], color=colors[0], label=r'$E_\mathrm{tot}$')
            axs_field.scatter(t[i], energy[0][i], color=colors[0])
            axs_field.set_xlim(t[0], t[-1])
            axs_field.set_xlabel(f'$t$ ({time_unit:s})')
            axs_field.set_ylabel(r'$E$ (a.u.)')

        plt.tight_layout()
        # for gif
        if gif:
            qa.gif.save_frame(gif_frames, dpi_gif=dpi_gif)

        plt.show()
        plt.pause(pause)


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
if (dynamics == 'rt') and 'field_coupling' in namelist['rt']:
    use_field = namelist['rt']['field_coupling']
else:
    use_field = False

# check if exact_factor and ef_gauge are in the input file
if 'exact_factor' in namelist['rt']:
    ef = namelist['rt']['exact_factor']
else:
    print("ERROR: Exact factorization keyword 'exact_factor' not found in the input file.")
    exit(1)

if 'ef_gauge' in namelist['rt']:
    ef_gauge = namelist['rt']['ef_gauge']
else:
    print("Default EF gauge considered.")
    ef_gauge = 'A0'



########## initialize ##########
if ef:
    print("\nPlotting exact factorization quantities from 1D quantum dynamics.\n")
else:
    print("ERROR: Exact factorization quantities not calculated! Exiting...")
    exit(1)

# reading wave function and all other important data
wf_bh, pot_en, [x, y, z], nframes_wf = qa.read.wf(rank=rank, nstates=nstates, xngrid=xngrid, yngrid=yngrid,
                                                  adiabatic=adiabatic)
# calculating Born-Huang densities from nuclear wave functions
den_bh = np.array(np.real(np.conjugate(wf_bh) * wf_bh), dtype=float)

# reading energies and time
en_states, nframes = qa.read.energies(dynamics=dynamics, nstates=nstates)
energy = en_states[0][1:]
t = en_states[0][0]

# reading populations
pop_ad, pop_diab = qa.read.pop()

# reading field
if use_field:
    field = qa.read.field()

# check consistency of the number of frames
if nframes != nframes_wf:
    print("Number of data in wf files and energy files is not matching. Exiting..")
    exit(1)

# reading exact factorization quantities
nucdens, nucphase, gitdpes, gdtdpes, tdvp, el_coeff, [ef_x, ef_y] = qa.read.ef(rank=rank, nstates=nstates,
                                                                               xngrid=xngrid, yngrid=yngrid)

# check consistency of the number of frames
if nframes != len(nucdens):
    print("Number of data in wf files and ef files is not matching. Exiting..")
    exit(1)

########## code ##########


########## plotting ##########
plt.rcParams["font.family"] = 'Helvetica'
# setting up online plotting
if plot_online:
    print("*plotting on-the-fly")
    plt.ion()
else:
    gif = False

# gifs for plotting
if gif:
    gif_frames = qa.gif.init_gif()

# converting time units
if time_unit == 'fs':
    pop_ad[0] *= qa.units.autofs
    pop_diab[0] *= qa.units.autofs
    en_states[:, 0] *= qa.units.autofs
    if use_field: field[0] *= qa.units.autofs

# scaling wf and density for plotting purposes
# note that this might cause troubles during analysis in not renormalized
wf_bh *= wf_scaling
den_bh *= wf_scaling ** 2
nucdens *= wf_scaling ** 2

if rank == 1:
    plot_ef_1d()

# if plot_on_the_fly:
if plot_online:
    if not gif:
        plt.pause(200.0)
    else:
        plt.pause(0.1)
    plt.ioff()

if gif:
    qa.gif.make_gif(gif_frames, duration, gif_name=f'ef_{ef_gauge}_{rank}D_{nstates}states')
