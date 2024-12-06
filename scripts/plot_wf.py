"""Plotting wave functions produced by QDyn. Versions for real-time and imaginary-time are available. Only 2D and 1D
wave functions can be plotted."""

from os.path import exists

import matplotlib.pyplot as plt
import numpy as np

import qdyn_analyze as qa

########## INPUT ##########
# general plotting settings
plot_online = True  # if true, plots the results on the fly, otherwise it's just the final plot
pause = 0.00001  # pause between frames during plotting
frame_step = 1  # plot only frame_step instead of plotting every frame
wf_scaling = 0.005  # scaling factor for the wf and density so that they are visible in the plots

# units for plotting
time_unit = ['a.u.', 'fs'][1]  # time units in femtoseconds or atomic time units

# real-time multistate options
adiabatic = True  # for rt dynamics with states more than 2

# norm of wave function
plot_norm = True
plot_wf = False

# 2D plots
heatmap = True

# gif
gif = False
dpi_gif = 100
duration = 2  # in ms, (duration=1000 * 1/fps, fps=frames per second).


########## functions ##########
def plot_rt_1d():
    # setting canvas
    if nstates == 1:
        fig, axs = plt.subplots(3, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [4, 1, 1]})
        axs_wf = axs[0]
        axs_en = axs[1]
        axs_field = axs[2]
    else:
        fig, axs = plt.subplots(4, 1, figsize=(5, 10), gridspec_kw={'height_ratios': [4, 1, 1, 1]})
        axs_wf = axs[0]
        axs_pop = axs[1]
        axs_en = axs[2]
        axs_field = axs[3]

    # getting ratios
    if plot_wf:
        vmax = 1.01 * np.max(wf)
        vmin = np.min([np.min(wf), np.min(pot_en)])
    else:
        vmax = 1.01 * np.max(den)
        vmin = np.min([np.min(den), np.min(pot_en)])
    emin = np.min(energy[:-2])
    emax = np.max(energy[:-2])
    emin -= 0.1 * (emax - emin)
    emax += 0.1 * (emax - emin)

    if plot_online:
        for i in range(0, nframes, frame_step):
            axs_en.cla()
            axs_wf.cla()
            for j in range(0, nstates):

                axs_wf.plot(x, pot_en[j], color='black')

                if plot_norm:
                    axs_wf.plot(x, den[j][i] + pot_en[j], linewidth=0.5, label=rf'$|\psi_{j:d}|^2$')
                    axs_wf.fill_between(x, den[j][i] * 0 + pot_en[j], den[j][i] + pot_en[j], alpha=0.35)
                if plot_wf:
                    axs_wf.plot(x, np.real(wf[j][i]) + pot_en[j], linewidth=0.5, label=rf'Re[$\psi_{j:d}$]')
                    axs_wf.plot(x, np.imag(wf[j][i]) + pot_en[j], linewidth=0.5, label=rf'Im[$\psi_{j:d}$]')
                    if not plot_norm:
                        axs_wf.fill_between(x, np.real(wf[j][i]) * 0 + pot_en[j], np.real(wf[j][i]) + pot_en[j],
                                            alpha=0.15)
                        axs_wf.fill_between(x, np.imag(wf[j][i]) * 0 + pot_en[j], np.imag(wf[j][i]) + pot_en[j],
                                            alpha=0.15)

            axs_wf.set_xlim(xmin, xmax)
            axs_wf.set_ylim(vmin, vmax)
            axs_wf.set_ylabel(r'$E$ (a.u.)')
            axs_wf.set_xlabel(r'$x$ (a.u.)')
            axs_wf.legend(labelspacing=0)

            axs_en.plot(t[:i + 1], energy[0][:i + 1], color='black', label=r'$E_\mathrm{tot}$')
            axs_en.scatter(t[i], energy[0][i], color='black')
            axs_en.plot(t[:i + 1], energy[1][:i + 1], color='black', linestyle='dashed', label=r'$V$')
            axs_en.scatter(t[i], energy[1][i], color='black')
            axs_en.plot(t[:i + 1], energy[2][:i + 1], color='black', linestyle='dotted', label=r'$E_\mathrm{k}$')
            axs_en.scatter(t[i], energy[2][i], color='black')

            axs_en.set_xlim(t[0], t[-1])
            axs_en.set_ylim(emin, emax)
            axs_en.set_ylabel(r'$E$ (a.u.)')
            axs_en.legend()

            axs_wf.set_title(f'time: {t[i]:.0f} {time_unit:s}  $E={energy[0, i]:.4f}$ a.u.')
            axs_en.set_xlabel(f'$t$ ({time_unit:s})')

            if use_field:
                axs_field.cla()
                axs_field.plot(field[0, :i], field[1, :i])
                axs_field.set_xlim(field[0, 0], field[0, -1])
                axs_field.set_xlim(t[0], t[-1])
                axs_field.set_xlabel(f'$t$ ({time_unit:s})')

            else:
                axs_field.cla()
                axs_field.plot(t[:i + 1], energy[0][:i + 1], color='black', label=r'$E_\mathrm{tot}$')
                axs_field.scatter(t[i], energy[0][i], color='black')
                axs_field.set_xlim(t[0], t[-1])
                axs_field.set_xlabel(f'$t$ ({time_unit:s})')
                axs_field.set_ylabel(r'$E$ (a.u.)')

            if nstates > 1:
                axs_pop.cla()
                for state in range(1, nstates + 1):
                    if adiabatic:
                        axs_pop.plot(pop_ad[0, :i], pop_ad[state, :i], label=f'state {state:d}')
                        axs_pop.set_ylabel('Adiabitc populaitons')
                    else:
                        axs_pop.plot(pop_diab[0], pop_diab[state], label=f'state {state:d}')
                        axs_pop.set_ylabel('Diabatic populaitons')
                if adiabatic:
                    axs_pop.plot(pop_ad[0, :i], pop_ad[-1, :i], label=f'norm', color='black')
                else:
                    axs_pop.plot(pop_diab[0, :i], pop_diab[-1, :i], label=f'norm', color='black')

                axs_pop.set_xlim(t[0], t[-1])
                axs_pop.set_xlabel(f'$t$ ({time_unit:s})')
                axs_pop.legend(frameon=False, labelspacing=0)

            plt.tight_layout()
            # for gif
            if gif:
                qa.gif.save_frame(gif_frames, dpi_gif=dpi_gif)

            if plot_online:
                plt.show()
                plt.pause(pause)
    else:
        for j in range(0, nstates):
            axs_wf.plot(x, pot_en[j], color='black')

            if plot_norm:
                axs_wf.plot(x, den[j][-1] + pot_en[j], linewidth=0.5, label=rf'$|\psi_{j:d}|^2$')
                axs_wf.fill_between(x, den[j][-1] * 0 + pot_en[j], den[j][-1] + pot_en[j], alpha=0.35)
            if plot_wf:
                axs_wf.plot(x, np.real(wf[j][-1]) + pot_en[j], linewidth=0.5, label=rf'Re[$\psi_{j:d}$]')
                axs_wf.plot(x, np.imag(wf[j][-1]) + pot_en[j], linewidth=0.5, label=rf'Im[$\psi_{j:d}$]')
                if not plot_norm:
                    axs_wf.fill_between(x, np.real(wf[j][-1]) * 0 + pot_en[j], np.real(wf[j][-1]) + pot_en[j],
                                        alpha=0.15)
                    axs_wf.fill_between(x, np.imag(wf[j][-1]) * 0 + pot_en[j], np.imag(wf[j][-1]) + pot_en[j],
                                        alpha=0.15)

        axs_wf.set_xlim(xmin, xmax)
        axs_wf.set_ylim(vmin, vmax)
        axs_wf.set_ylabel(r'$E$ (a.u.)')
        axs_wf.set_xlabel(r'$x$ (a.u.)')
        axs_wf.legend(labelspacing=0)

        axs_en.plot(t, energy[0], color='black', label=r'$E_\mathrm{tot}$')
        axs_en.plot(t, energy[1], color='black', linestyle='dashed', label=r'$V$')
        axs_en.plot(t, energy[2], color='black', linestyle='dotted', label=r'$E_\mathrm{k}$')

        axs_en.set_xlim(t[0], t[-1])
        axs_en.set_ylim(emin, emax)
        axs_en.set_ylabel(r'$E$ (a.u.)')
        axs_en.legend()

        axs_wf.set_title(f'time: {t[-1]:.0f} {time_unit:s}  $E={energy[0, -1]:.4f}$ a.u.')
        axs_en.set_xlabel(f'$t$ ({time_unit:s})')

        if use_field:
            axs_field.plot(field[0], field[1])
            axs_field.set_xlim(field[0, 0], field[0, -1])
            axs_field.set_xlim(t[0], t[-1])
            axs_field.set_xlabel(f'$t$ ({time_unit:s})')
            axs_field.set_ylabel(r'$\vec{E}(t)$ / a.u.')
        else:
            axs_field.plot(t, energy[0], color='black', label=r'$E_\mathrm{tot}$')
            axs_field.set_xlim(t[0], t[-1])
            axs_field.set_xlabel(f'$t$ ({time_unit:s})')
            axs_field.set_ylabel(r'$E$ (a.u.)')

        if nstates > 1:
            for state in range(1, nstates + 1):
                if adiabatic:
                    axs_pop.plot(pop_ad[0], pop_ad[state], label=f'state {state:d}')
                    axs_pop.set_ylabel('Adiabitc populaitons')
                else:
                    axs_pop.plot(pop_diab[0], pop_diab[state], label=f'state {state:d}')
                    axs_pop.set_ylabel('Diabatic populaitons')
            if adiabatic:
                axs_pop.plot(pop_ad[0], pop_ad[-1], label=f'norm', color='black')
            else:
                axs_pop.plot(pop_diab[0], pop_diab[-1], label=f'norm', color='black')

            axs_pop.set_xlim(t[0], t[-1])
            axs_pop.set_xlabel(f'$t$ ({time_unit:s})')
            axs_pop.legend(frameon=False, labelspacing=0)

        plt.tight_layout()
        plt.show()


def plot_it_1d():
    # setting canvas
    colors = plt.cm.viridis(np.linspace(0, 1, nstates))
    plt.rcParams["font.family"] = 'Helvetica'

    # setting up figures
    if final_plot:
        fig, axs = plt.subplots(1, 1, figsize=(5, 6.5))
        axs.minorticks_on()
        axs.tick_params('both', direction='in', which='both', top=True, right=True)
        axs.tick_params('y', direction='in', which='minor', length=0)
    else:
        fig, axs = plt.subplots(2, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [2, 1]})

    # setting up maximum and minimum values of axes for the two different scenarios
    if plot_wf:
        vmax = 1.0 * np.max([np.max(np.real(wf[-1, :])), np.max(np.imag(wf[-1, :]))]) + np.max(en_states[:, 1, -1])
        vmin = np.min([np.max(np.real(wf[-1, :])) + np.min(en_states[:, 1, -1]),
                       np.max(np.imag(wf[-1, :])) + np.min(en_states[:, 1, -1]), np.min(pot_en)])
        sigma_e = np.std(en_states[:, 1, :]) * 5
        aver_e = np.average(en_states[:, 1:3, :])
        emin = np.min([np.min(en_states[:, 1:3, :]), aver_e - sigma_e])
        emax = np.max([np.max(en_states[:, 1:3, :]), aver_e + sigma_e])
        print("vmin,vmax,emin,emax:", vmin, vmax, emin, emax)
    else:  # when plotting norm
        vmax = 1.5 * np.max(den[-1, -1]) + en_states[-1, 1, -1]
        vmin = np.min(pot_en)
        sigma_e = np.std(en_states[:, 1, :]) * 5
        aver_e = np.average(en_states[:, 1:3, :])
        emin = np.min([np.min(en_states[:, 1:3, :]), aver_e - sigma_e])
        emax = np.max([np.max(en_states[:, 1:3, :]), aver_e + sigma_e])
        print("vmin,vmax,emin,emax:", vmin, vmax, emin, emax)

    if plot_online:
        for j in range(0, nstates):  # loop over all states
            for i in range(1, nframes, frame_step):  # going through frames
                axs[0].cla()
                axs[1].cla()

                axs[0].plot(x, pot_en[j], color='black')

                # plotting previous states
                for k in range(0, j):
                    if plot_norm:
                        axs[0].plot(x, den[k, -1] + en_states[k, 1, -1], linewidth=0.5, color=colors[k])
                        axs[0].fill_between(x, x * 0 + en_states[k, 1, -1], den[k, -1] + en_states[k, 1, -1],
                                            alpha=0.15, color=colors[k])
                    elif plot_wf:
                        axs[0].fill_between(x, x * 0 + en_states[k, 1, -1], np.real(wf[k, -1]) + en_states[k, 1, -1],
                                            alpha=0.15, color=colors[k])
                        axs[0].fill_between(x, x * 0 + en_states[k, 1, -1], np.imag(wf[k, -1]) + en_states[k, 1, -1],
                                            alpha=0.15, color=colors[nstates - k - 1])

                    # plot state index on the right of the density/wf
                    axs[0].text(7.1, en_states[k, 1, -1], r'$n=\mathregular{%d}$' % (k + 1), horizontalalignment='left')

                    # plot total energy in the bottom plot
                    axs[1].plot(en_states[k, 0, :], en_states[k, 1, :], alpha=1, color=colors[k])

                if plot_norm:
                    axs[0].plot(x, den[j, i] + en_states[j, 1, i], linewidth=0.5, label=r'$|\psi_{%d}|^2$' % j,
                                color=colors[j])
                    axs[0].fill_between(x, den[j, i] * 0 + en_states[j, 1, i], den[j, i] + en_states[j, 1, i],
                                        alpha=0.35, color=colors[j])
                if plot_wf:
                    axs[0].plot(x, np.real(wf[j, i]) + en_states[j, 1, i], linewidth=0.5, label=r'Re{$\psi$}',
                                color=colors[j])
                    axs[0].plot(x, np.imag(wf[j, i]) + en_states[j, 1, i], linewidth=0.5, label=r'Im{$\psi$}',
                                color=colors[nstates - j - 1])
                    if not plot_norm:
                        axs[0].fill_between(x, np.real(wf[j, i]) * 0 + en_states[j, 1, i],
                                            np.real(wf[j, i]) + en_states[j, 1, i], alpha=0.15, color=colors[j])
                        axs[0].fill_between(x, np.imag(wf[j, i]) * 0 + en_states[j, 1, i],
                                            np.imag(wf[j, i]) + en_states[j, 1, i], alpha=0.15,
                                            color=colors[nstates - j - 1])

                # plot state index on the right of the density/wf
                axs[0].text(7.1, en_states[j, 1, -1], r'$n=\mathregular{%d}$' % (j + 1), horizontalalignment='left')

                axs[0].set_ylim(vmin, vmax)
                axs[0].set_ylabel(r'$E$ (a.u.)')
                axs[0].set_xlabel(r'$x$ (a.u.)')
                axs[0].legend(labelspacing=0)

                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 1, :i + 1], color='black', label=r'$E_\mathrm{tot}$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 1, i], color='black')
                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 2, :i + 1], color='black', linestyle='dashed',
                            label=r'$V$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 2, i], color='black')
                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 3, :i + 1], color='black', linestyle='dotted',
                            label=r'$E_\mathrm{k}$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 3, i], color='black')
                axs[1].set_xlim(en_states[j, 0, 0], en_states[j, 0, -1])
                axs[1].set_ylim(emin, emax)
                axs[1].set_ylabel(r'$E$ (a.u.)')
                axs[1].legend()

                current_time = en_states[j, 0, i]
                current_energy = en_states[j, 1, i]
                axs[0].set_title(f'$t={current_time:.0f}$ {time_unit:s};  $E={current_energy:.8f}$ a.u.')
                axs[1].set_xlabel(f'$t$ ({time_unit:s})')

                plt.tight_layout()
                # for gif
                if gif:
                    qa.gif.save_frame(gif_frames, dpi_gif=dpi_gif)

                # end gif
                plt.show()
                plt.pause(pause)
    else:  # final plot
        for j in range(nstates):
            axs.plot(x, pot_en[j], color='black')
            axs.text(7.1, en_states[j, 1, -1], r'$n=\mathregular{%d}$' % (j + 1), horizontalalignment='left')

            if plot_norm:
                axs.plot(x, den[j, -1] + en_states[j, 1][-1], linewidth=0.5, label=r'$|\psi|^2$', color=colors[j])
                axs.fill_between(x, den[j, -1] * 0 + en_states[j, 1][-1], den[j, -1] + en_states[j, 1][-1], alpha=0.35,
                                 color=colors[j])
            if plot_wf:
                axs.plot(x, np.real(wf[j, -1]) + en_states[j, 1][-1], linewidth=0.5, label=r'Re{$\psi$}',
                         color=colors[j])
                axs.plot(x, np.imag(wf[j, -1]) + en_states[j, 1][-1], linewidth=0.5, label=r'Im{$\psi$}',
                         color=colors[nstates - j - 1])
                if not plot_norm:
                    axs.fill_between(x, np.real(wf[j, -1]) * 0 + en_states[j, 1][-1],
                                     np.real(wf[j, -1]) + en_states[j, 1][-1], alpha=0.15, color=colors[j])
                    axs.fill_between(x, np.imag(wf[j, -1]) * 0 + en_states[j, 1][-1],
                                     np.imag(wf[j, -1]) + en_states[j, 1][-1], alpha=0.15,
                                     color=colors[nstates - j - 1])

        axs.set_ylim(bottom=vmin, top=vmax)
        axs.set_xlim(left=xmin, right=xmax)
        axs.set_ylabel(r'$E$ (a.u.)')
        axs.set_xlabel(r'$x$ (a.u.)')
        axs.set_yticks(en_states[:, 1, -1])
        plt.tight_layout()
        plt.savefig("it_1d_eigenstates", dpi=300)
        plt.show()


# todo: modify
def plot_it_2d():
    # setting canvas
    if final_plot:
        fig = plt.figure(figsize=(5, 5))
        axs = [fig.add_subplot(111, projection='3d')]
    else:
        fig = plt.figure(figsize=(5, 7.5))
        axs = [fig.add_subplot(211, projection='3d'), fig.add_subplot(212)]

    if final_plot:
        for j in range(nstates):
            axs[0].plot_trisurf(wf[j, -1, 0, :], wf[j, -1, 1, :], wf[j, -1, 2, :])
        axs[0].set_xlabel('$x$ (a.u.)')
        axs[0].set_ylabel('y / a.u.')
        axs[0].set_zlabel('z / a.u.')
        plt.tight_layout()
        plt.show()
    else:
        if plot_wf:
            sigma_e = np.std(en_states[:, 1, :]) * 5
            aver_e = np.average(en_states[:, 1:3, :])
            emin = np.min([np.min(en_states[:, 1:3, :]), aver_e - sigma_e])
            emax = np.max([np.max(en_states[:, 1:3, :]), aver_e + sigma_e])
        else:  # when plotting norm
            sigma_e = np.std(en_states[:, 1, :]) * 5
            aver_e = np.average(en_states[:, 1:3, :])
            emin = np.min([np.min(en_states[:, 1:3, :]), aver_e - sigma_e])
            emax = np.max([np.max(en_states[:, 1:3, :]), aver_e + sigma_e])

        for j in range(0, nstates):
            for i in range(1, nframes):
                axs[0].cla()
                axs[1].cla()

                axs[0].plot_trisurf(wf[j, i, 0, :], wf[j, i, 1, :], wf[j, i, 2, :])
                # axs[0].plot_trisurf(wf[j, i, 0, :], wf[j, i, 1, :], wf[j, i, 3, :], alpha=0.5)

                axs[0].set_zlabel(r'$E$ (a.u.)')
                axs[0].set_ylabel(r'y / a.u.')
                axs[0].set_xlabel(r'$x$ (a.u.)')

                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 1, :i + 1], color='black', label=r'$E_\mathrm{tot}$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 1, i], color='black')
                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 2, :i + 1], color='black', linestyle='dashed',
                            label=r'$V$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 2, i], color='black')
                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 3, :i + 1], color='black', linestyle='dotted',
                            label=r'$E_\mathrm{k}$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 3, i], color='black')
                axs[1].set_xlim(en_states[j, 0, 0], en_states[j, 0, -1])
                axs[1].set_ylim(emin, emax)
                axs[1].set_ylabel(r'$E$ (a.u.)')
                axs[1].legend()

                axs[0].set_title('imaginary time: %d a.u.  $E=%.4f$a.u.' % (en_states[j, 0, i], en_states[j, 1, i]))
                axs[1].set_xlabel(r'$t$ (a.u.)')

                plt.tight_layout()
                # for gif
                if gif:
                    qa.gif.save_frame(gif_frames, dpi_gif=dpi_gif)

                # end gif
                if plot_online:
                    plt.show()
                    plt.pause(pause)


# todo: modify
def plot_it_2d_heatmap():
    x, y = np.meshgrid(wf[0, -1, 0, ::xngrid], wf[0, -1, 1, :yngrid])

    if final_plot:
        for j in range(nstates):
            fig, axs = plt.subplots(1, 1, figsize=(6, 5))

            if plot_norm:
                z = np.reshape(wf[j, -1, 4, :], (xngrid, yngrid)).T
                pc = axs.pcolormesh(x, y, z, cmap='viridis')
            else:
                z = np.reshape(wf[j, -1, 2, :], (xngrid, yngrid)).T
                if j == 0:
                    pc = axs.pcolormesh(x, y, z, cmap='Blues')
                else:
                    pc = axs.pcolormesh(x, y, z, cmap='RdBu')

            fig.colorbar(pc, ax=axs)

            axs.set_ylabel(r'y / a.u.')
            axs.set_xlabel(r'$x$ (a.u.)')
            axs.set_title(r'$E=%.3f a.u.$' % en_states[j, 1, -1])

            plt.tight_layout()
            plt.show()
    else:

        fig, axs = plt.subplots(2, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [2, 1]})

        emin = np.min(en_states[:, 1:3, :])
        emax = np.max(en_states[:, 1:3, :])
        emin -= 0.1 * (emax - emin)
        emax += 0.1 * (emax - emin)

        for j in range(0, nstates):
            for i in range(1, nframes):
                axs[0].cla()
                axs[1].cla()

                if plot_norm:
                    z = np.reshape(wf[j, i, 4, :], (xngrid, yngrid)).T
                    axs[0].pcolormesh(x, y, z, cmap='viridis')
                else:
                    z = np.reshape(wf[j, i, 2, :], (xngrid, yngrid)).T
                    axs[0].pcolormesh(x, y, z, cmap='RdBu')

                axs[0].set_ylabel(r'y / a.u.')
                axs[0].set_xlabel(r'$x$ (a.u.)')

                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 1, :i + 1], color='black', label=r'$E_\mathrm{tot}$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 1, i], color='black')
                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 2, :i + 1], color='black', linestyle='dashed',
                            label=r'$V$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 2, i], color='black')
                axs[1].plot(en_states[j, 0, :i + 1], en_states[j, 3, :i + 1], color='black', linestyle='dotted',
                            label=r'$E_\mathrm{k}$')
                axs[1].scatter(en_states[j, 0, i], en_states[j, 3, i], color='black')
                axs[1].set_xlim(en_states[j, 0, 0], en_states[j, 0, -1])
                axs[1].set_ylim(emin, emax)
                axs[1].set_ylabel(r'$E$ (a.u.)')
                axs[1].legend()

                axs[0].set_title('imaginary time: %d a.u.  $E=%.4f$a.u.' % (en_states[j, 0, i], en_states[j, 1, i]))
                axs[1].set_xlabel(r'$t$ (a.u.)')

                plt.tight_layout()
                # for gif
                if gif:
                    qa.gif.save_frame(gif_frames, dpi_gif=dpi_gif)

                # end gif
                if plot_online:
                    plt.show()
                    plt.pause(pause)


# todo: modify
def plot_rt_2d_heatmap():
    fig, axs = plt.subplots(2, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [2, 1]})

    emin = np.min(energy[:-2])
    emax = np.max(energy[:-2])
    emin -= 0.1 * (emax - emin)
    emax += 0.1 * (emax - emin)

    x, y = np.meshgrid(wf[0, -1, 0, ::xngrid], wf[0, -1, 1, :yngrid])

    for j in range(0, nstates):
        for i in range(1, nframes):
            axs[0].cla()
            axs[1].cla()

            if plot_norm:
                z = np.reshape(wf[j, i, 4, :], (xngrid, yngrid))
                axs[0].pcolormesh(x, y, z, cmap='viridis')
            else:
                z = np.reshape(wf[j, i, 2, :], (xngrid, yngrid))
                axs[0].pcolormesh(x, y, z, cmap='RdBu')

            axs[0].set_ylabel(r'y / a.u.')
            axs[0].set_xlabel(r'$x$ (a.u.)')

            axs[1].plot(t[:i + 1], energy[0][:i + 1], color='black', label=r'$E_\mathrm{tot}$')
            axs[1].scatter(t[i], energy[0][i], color='black')
            axs[1].plot(t[:i + 1], energy[1][:i + 1], color='black', linestyle='dashed', label=r'$V$')
            axs[1].scatter(t[i], energy[1][i], color='black')
            axs[1].plot(t[:i + 1], energy[2][:i + 1], color='black', linestyle='dotted', label=r'$E_\mathrm{k}$')
            axs[1].scatter(t[i], energy[2][i], color='black')
            axs[1].set_xlim(t[0], t[-1])
            axs[1].set_ylim(emin, emax)
            axs[1].set_ylabel(r'$E$ (a.u.)')
            axs[1].legend()

            axs[0].set_title('time: %d a.u.  $E=%.4f$a.u.' % (t[i], energy[0, i]))
            axs[1].set_xlabel(r'$t$ (a.u.)')

            plt.tight_layout()
            # for gif
            if gif:
                qa.gif.save_frame(gif_frames, dpi_gif=dpi_gif)

            # end gif
            if plot_online:
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

########## initialize ##########
if dynamics == 'rt':
    print("\nPlotting real time dynamics from Qdyn.\n")
    if nstates == 1:
        adiabatic = False  # only diabatic wf ploted because it is the same as adiabatic
elif dynamics == 'it':
    print("\nPlotting imaginary time dynamics from Qdyn.\n")
    adiabatic = False

if rank == 3:
    print("ERROR: 3D plotting is not available! Exiting...")
    exit(1)

# reading wave function and all other important data
wf, pot_en, [x, y, z], nframes_wf = qa.read.wf(rank=rank, nstates=nstates, xngrid=xngrid, yngrid=yngrid,
                                               adiabatic=adiabatic)
en_states, nframes = qa.read.energies(dynamics=dynamics, nstates=nstates)

# calculating density from wave function
den = np.array(np.real(np.conjugate(wf) * wf), dtype=float)

if use_field:
    field = qa.read.field()

if time_unit == 'fs':
    en_states[:, 0] *= qa.units.autofs
    if use_field: field[0] *= qa.units.autofs

if dynamics == 'rt':
    energy = en_states[0][1:]
    t = en_states[0][0]
    if nstates > 1:
        pop_ad, pop_diab = qa.read.pop()
        if time_unit == 'fs':
            pop_ad[0] *= qa.units.autofs
            pop_diab[0] *= qa.units.autofs

# estimate number of frames
if nframes != nframes_wf:
    print("Number of data in wf files and energy files is not matching. Exiting..")
    exit(1)

# todo: this must be removed once plot_final is completely gone
final_plot = False
if not plot_online:
    final_plot = True

if plot_online:
    print("*plotting on-the-fly")
    plt.ion()
else:
    gif = False

if gif:
    gif_frames = qa.gif.init_gif()

########## code ##########

# scaling wf and density for plotting purposes
# note that this might cause troubles during analysis in not renormalized
wf *= wf_scaling
den *= wf_scaling

if rank == 1:
    # plot wf
    if dynamics == 'rt':
        plot_rt_1d()
    elif dynamics == 'it':
        plot_it_1d()
elif rank == 2:
    if dynamics == 'rt':
        if heatmap:
            plot_rt_2d_heatmap()
    elif dynamics == 'it':
        if heatmap:
            plot_it_2d_heatmap()
        else:
            plot_it_2d()

# if plot_on_the_fly:
if plot_online:
    plt.pause(50.0)
    plt.ioff()

if gif:
    qa.gif.make_gif(gif_frames, duration)
