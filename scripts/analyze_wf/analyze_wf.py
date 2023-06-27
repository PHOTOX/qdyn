# analyze wf from QDyn

from os import makedirs, remove, rmdir
from os.path import exists

import f90nml
import imageio
import matplotlib.pyplot as plt
import numpy as np

########## INPUT ##########
# plotting
final_plot = True
plot_online = True
plot_norm = False
plot_wf = False
heatmap = True
pause = 0.0001
# gif
gif = False
dpi_gif = 100
duration = 2  # in ms, (duration=1000 * 1/fps, fps=frames per second).


########## functions ##########
def read_wf():
    # TODO: add loop over states
    # TODO: Also add reading the file
    wf = []
    for i in range(0, nstates):
        if rank == 1:
            wf_file = 'wf1d.%d.out' % (i + 1)
            grid_size = ngrid
        elif rank == 2:
            wf_file = 'wf2d.%d.out' % (i + 1)
            grid_size = ngrid ** 2

        if exists(wf_file):
            wf.append(np.genfromtxt(wf_file))
            # estimate number of frames
            if np.shape(wf[i])[0] % grid_size != 0:
                print("Number of lines and number of grid point are not matching. Exiting..")
                exit(1)
            nframes_wf = int(np.shape(wf[i])[0] / grid_size)
            if rank == 1:
                wf[i] = np.reshape(wf[i], (nframes_wf, grid_size, 5)).transpose((0, 2, 1))
            elif rank == 2:
                wf[i] = np.reshape(wf[i], (nframes_wf, grid_size, 6)).transpose((0, 2, 1))
        else:
            print("File '%s' does not exist. Exiting.." % wf_file)
            exit(1)
    return np.array(wf), nframes_wf


def read_energies():
    en_states = []
    if run == 0:
        en_file = 'energies.dat'
        if exists(en_file):
            en_states.append(np.genfromtxt(en_file).T)
            print("'%s' read" % en_file)
        else:
            print("File '%s' does not exist. Exiting.." % en_file)
            exit(1)
        nframes = np.shape(en_states)[2]
    elif run == 1:
        for i in range(0, nstates):
            en_file = 'energies.%d.dat' % (i + 1)
            if exists(en_file):
                en_states.append(np.genfromtxt(en_file).T)
                print("'%s' read" % en_file)
            else:
                print("File '%s' does not exist. Exiting.." % en_file)
                exit(1)
        en_states = np.array(en_states)
        print("\nTotal energies:", en_states[:, 1, -1])
        nframes = np.shape(en_states)[2]
    return en_states, nframes


def read_field():
    return np.genfromtxt('field.dat').T


def gif_frame(i):
    fig_file = fig_folder + '/frame_%d.png' % i
    plt.savefig(fig_file, format='png', dpi=dpi_gif)
    frames.append(imageio.v2.imread(fig_file))
    remove(fig_file)


def plot_rt_1d():
    # setting canvas
    if use_field:
        fig, axs = plt.subplots(3, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [4, 1, 1]})
    else:
        fig, axs = plt.subplots(2, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [2, 1]})

    # getting ratios
    vmax = np.max([1.2 * np.max(wf[:, :, 1:3, :]), np.max(energy)])
    vmin = np.min(wf[:, :, 1:, :])
    emin = np.min(energy[:-2])
    emax = np.max(energy[:-2])
    emin -= 0.1 * (emax - emin)
    emax += 0.1 * (emax - emin)

    for j in range(0, nstates):
        for i in range(nframes):
            axs[0].cla()

            axs[0].plot(wf[j][0][0], wf[j][0][4], color='black')

            if plot_norm:
                axs[0].plot(wf[j][i][0], wf[j][i][3] + energy[0][i], linewidth=0.5, label=r'$|\psi|^2$')
                axs[0].fill_between(wf[j][i][0], wf[j][i][3] * 0 + energy[0][i], wf[j][i][3] + energy[0][i], alpha=0.35)
            if plot_wf:
                axs[0].plot(wf[j][i][0], wf[j][i][1] + energy[0][i], linewidth=0.5, label=r'Re{$\psi$}')
                axs[0].plot(wf[j][i][0], wf[j][i][2] + energy[0][i], linewidth=0.5, label=r'Im{$\psi$}')
                if not plot_norm:
                    axs[0].fill_between(wf[j][i][0], wf[j][i][1] * 0 + energy[0][i], wf[j][i][1] + energy[0][i],
                                        alpha=0.15)
                    axs[0].fill_between(wf[j][i][0], wf[j][i][2] * 0 + energy[0][i], wf[j][i][2] + energy[0][i],
                                        alpha=0.15)

            axs[0].set_ylim(vmin, vmax)
            axs[0].set_ylabel(r'E / a.u.')
            axs[0].set_xlabel(r'x / a.u.')
            axs[0].legend(labelspacing=0)

            axs[1].cla()

            axs[1].plot(t[:i + 1], energy[0][:i + 1], color='black', label=r'$E_\mathrm{tot}$')
            axs[1].scatter(t[i], energy[0][i], color='black')
            axs[1].plot(t[:i + 1], energy[1][:i + 1], color='black', linestyle='dashed', label=r'$V$')
            axs[1].scatter(t[i], energy[1][i], color='black')
            axs[1].plot(t[:i + 1], energy[2][:i + 1], color='black', linestyle='dotted', label=r'$E_\mathrm{k}$')
            axs[1].scatter(t[i], energy[2][i], color='black')
            axs[1].set_xlim(t[0], t[-1])
            # axs[1].set_ylim(emin, emax)
            axs[1].set_ylabel(r'E / a.u.')
            axs[1].legend()

            axs[0].set_title('time: %d a.u.  $E=%.4f$a.u.' % (t[i], energy[0,i]))
            axs[1].set_xlabel(r'$t$ / a.u.')

            if use_field:
                axs[2].cla()
                axs[2].plot(field[0, :i], field[1, :i])
                axs[2].set_xlim(field[0, 0], field[0, -1])

            plt.tight_layout()
            # for gif
            if gif:
                gif_frame(i)

            if plot_online:
                plt.show()
                plt.pause(pause)


def plot_it_1d():
    # setting canvas
    if final_plot:
        fig, axs = plt.subplots(1, 1, figsize=(5, 5))
    else:
        fig, axs = plt.subplots(2, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [2, 1]})

    # TODO: problem with max value
    if plot_wf:
        vmax = np.max([1.2 * np.max(wf[-1, :, 1:3, :]) + np.max(en_states[:, 1, :])])
        # print(wf[:,:, 1:, :])
        vmin = np.min([np.min(wf[:, :, 1:3, :]) + np.min(en_states[:, 1, :]), np.min(wf[:, :, 4, :])])
        # vmin = -1
        sigma_e = np.std(en_states[:, 1, :]) * 5
        aver_e = np.average(en_states[:, 1:3, :])
        emin = np.min([np.min(en_states[:, 1:3, :]), aver_e - sigma_e])
        emax = np.max([np.max(en_states[:, 1:3, :]), aver_e + sigma_e])
        print("vmin,vmax,emin,emax:", vmin, vmax, emin, emax)
    else:  # when plotting norm
        print('here')
        vmax = 1.5 * np.max(wf[-1, -1, 3, :]) + np.max(en_states[-1, 1, -1])
        # print(wf[:,:, 1:, :])
        vmin = np.min(wf[0, 0, 4, :]) - 0.05
        # vmin = -1
        sigma_e = np.std(en_states[:, 1, :]) * 5
        aver_e = np.average(en_states[:, 1:3, :])
        emin = np.min([np.min(en_states[:, 1:3, :]), aver_e - sigma_e])
        emax = np.max([np.max(en_states[:, 1:3, :]), aver_e + sigma_e])
        print("vmin,vmax,emin,emax:", vmin, vmax, emin, emax)

    if final_plot:
        for j in range(nstates):
            axs.plot(wf[j][0][0], wf[j][0][4], color='black')
            axs.text(wf[j, 0, 0, 0], en_states[j, 1, -1], r'$n=%d$' % j)
            axs.text(wf[j, 0, 0, -1], en_states[j, 1, -1], r'$E=%.3f$' % en_states[j, 1, -1],
                     horizontalalignment='right')

            if plot_norm:
                axs.plot(wf[j][-1][0], wf[j][-1][3] + en_states[j, 1][-1], linewidth=0.5, label=r'$|\psi|^2$')
                axs.fill_between(wf[j][-1][0], wf[j][-1][3] * 0 + en_states[j, 1][-1],
                                 wf[j][-1][3] + en_states[j, 1][-1],
                                 alpha=0.35)
            if plot_wf:
                axs.plot(wf[j][-1][0], wf[j][-1][1] + en_states[j, 1][-1], linewidth=0.5, label=r'Re{$\psi$}')
                axs.plot(wf[j][-1][0], wf[j][-1][2] + en_states[j, 1][-1], linewidth=0.5, label=r'Im{$\psi$}')
                if not plot_norm:
                    axs.fill_between(wf[j][-1][0], wf[j][-1][1] * 0 + en_states[j, 1][-1],
                                     wf[j][-1][1] + en_states[j, 1][-1], alpha=0.15)
                    axs.fill_between(wf[j][-1][0], wf[j][-1][2] * 0 + en_states[j, 1][-1],
                                     wf[j][-1][2] + en_states[j, 1][-1], alpha=0.15)

        axs.set_ylim(bottom=vmin, top=vmax)
        axs.set_ylabel(r'E / a.u.')
        axs.set_xlabel(r'x / a.u.')
        # axs.legend(labelspacing=0)
        plt.tight_layout()
        plt.show()
    else:
        for j in range(0, nstates):
            for i in range(1, nframes):
                axs[0].cla()
                axs[1].cla()

                axs[0].plot(wf[j, 0, 0], wf[j, 0, 4], color='black')

                # ploting previous states
                # TODO: expand the plot previous function
                previous = False
                if previous:
                    for k in range(0, j):
                        axs[0].plot(wf[k, -1, 0], wf[k, -1, 3] + en_states[k, 1, -1], linewidth=0.5,
                                    label=r'$|\psi|^2$')
                        axs[0].fill_between(wf[k, -1, 0], wf[k, -1, 3] * 0 + en_states[k, 1, -1],
                                            wf[k, -1, 3] + en_states[k, 1, -1], alpha=0.15)
                        axs[1].plot(en_states[k, 0, :], en_states[k, 1, :], alpha=0.3)

                if plot_norm:
                    axs[0].plot(wf[j, i, 0], wf[j, i, 3] + en_states[j, 1, i], linewidth=0.5, label=r'$|\psi|^2$')
                    axs[0].fill_between(wf[j, i, 0], wf[j, i, 3] * 0 + en_states[j, 1, i],
                                        wf[j, i, 3] + en_states[j, 1, i], alpha=0.35)
                if plot_wf:
                    axs[0].plot(wf[j, i, 0], wf[j, i, 1] + en_states[j, 1, i], linewidth=0.5, label=r'Re{$\psi$}')
                    axs[0].plot(wf[j, i, 0], wf[j, i, 2] + en_states[j, 1, i], linewidth=0.5, label=r'Im{$\psi$}')
                    if not plot_norm:
                        axs[0].fill_between(wf[j, i, 0], wf[j, i, 1] * 0 + en_states[j, 1, i],
                                            wf[j, i, 1] + en_states[j, 1, i], alpha=0.15)
                        axs[0].fill_between(wf[j, i, 0], wf[j, i, 2] * 0 + en_states[j, 1, i],
                                            wf[j, i, 2] + en_states[j, 1, i], alpha=0.15)

                axs[0].set_ylim(vmin, vmax)
                axs[0].set_ylabel(r'E / a.u.')
                axs[0].set_xlabel(r'x / a.u.')
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
                axs[1].set_ylabel(r'E / a.u.')
                axs[1].legend()

                axs[0].set_title('imaginary time: %d a.u.  $E=%.4f$a.u.' % (en_states[j, 0, i], en_states[j, 1, i]))
                axs[1].set_xlabel(r'$t$ / a.u.')

                plt.tight_layout()
                # for gif
                if gif:
                    gif_frame(i)

                # end gif
                if plot_online:
                    plt.show()
                    plt.pause(pause)


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
        axs[0].set_xlabel('x / a.u.')
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

                axs[0].set_zlabel(r'E / a.u.')
                axs[0].set_ylabel(r'y / a.u.')
                axs[0].set_xlabel(r'x / a.u.')

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
                axs[1].set_ylabel(r'E / a.u.')
                axs[1].legend()

                axs[0].set_title('imaginary time: %d a.u.  $E=%.4f$a.u.' % (en_states[j, 0, i], en_states[j, 1, i]))
                axs[1].set_xlabel(r'$t$ / a.u.')

                plt.tight_layout()
                # for gif
                if gif:
                    gif_frame(i)

                # end gif
                if plot_online:
                    plt.show()
                    plt.pause(pause)


def plot_it_2d_heatmap():
    x, y = np.meshgrid(wf[0, -1, 0, ::ngrid], wf[0, -1, 1, :ngrid])

    if final_plot:
        for j in range(nstates):
            fig, axs = plt.subplots(1, 1, figsize=(6, 5))

            if plot_norm:
                z = np.reshape(wf[j, -1, 4, :], (ngrid, ngrid)).T
                pc = axs.pcolormesh(x, y, z, cmap='viridis')
            else:
                z = np.reshape(wf[j, -1, 2, :], (ngrid, ngrid)).T
                if j == 0:
                    pc = axs.pcolormesh(x, y, z, cmap='Blues')
                else:
                    pc = axs.pcolormesh(x, y, z, cmap='RdBu')

            fig.colorbar(pc,ax=axs)

            axs.set_ylabel(r'y / a.u.')
            axs.set_xlabel(r'x / a.u.')
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
                    z = np.reshape(wf[j, i, 4, :], (ngrid, ngrid)).T
                    axs[0].pcolormesh(x, y, z, cmap='viridis')
                else:
                    z = np.reshape(wf[j, i, 2, :], (ngrid, ngrid)).T
                    axs[0].pcolormesh(x, y, z, cmap='RdBu')

                axs[0].set_ylabel(r'y / a.u.')
                axs[0].set_xlabel(r'x / a.u.')

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
                axs[1].set_ylabel(r'E / a.u.')
                axs[1].legend()

                axs[0].set_title('imaginary time: %d a.u.  $E=%.4f$a.u.' % (en_states[j, 0, i], en_states[j, 1, i]))
                axs[1].set_xlabel(r'$t$ / a.u.')

                plt.tight_layout()
                # for gif
                if gif:
                    gif_frame(i)

                # end gif
                if plot_online:
                    plt.show()
                    plt.pause(pause)


def plot_rt_2d_heatmap():
    fig, axs = plt.subplots(2, 1, figsize=(5, 7.5), gridspec_kw={'height_ratios': [2, 1]})

    emin = np.min(energy[:-2])
    emax = np.max(energy[:-2])
    emin -= 0.1 * (emax - emin)
    emax += 0.1 * (emax - emin)

    x, y = np.meshgrid(wf[0, -1, 0, ::ngrid], wf[0, -1, 1, :ngrid])

    for j in range(0, nstates):
        for i in range(1, nframes):
            axs[0].cla()
            axs[1].cla()

            if plot_norm:
                z = np.reshape(wf[j, i, 4, :], (ngrid, ngrid))
                axs[0].pcolormesh(x, y, z, cmap='viridis')
            else:
                z = np.reshape(wf[j, i, 2, :], (ngrid, ngrid))
                axs[0].pcolormesh(x, y, z, cmap='RdBu')

            axs[0].set_ylabel(r'y / a.u.')
            axs[0].set_xlabel(r'x / a.u.')

            axs[1].plot(t[:i + 1], energy[0][:i + 1], color='black', label=r'$E_\mathrm{tot}$')
            axs[1].scatter(t[i], energy[0][i], color='black')
            axs[1].plot(t[:i + 1], energy[1][:i + 1], color='black', linestyle='dashed', label=r'$V$')
            axs[1].scatter(t[i], energy[1][i], color='black')
            axs[1].plot(t[:i + 1], energy[2][:i + 1], color='black', linestyle='dotted', label=r'$E_\mathrm{k}$')
            axs[1].scatter(t[i], energy[2][i], color='black')
            axs[1].set_xlim(t[0], t[-1])
            axs[1].set_ylim(emin, emax)
            axs[1].set_ylabel(r'E / a.u.')
            axs[1].legend()

            axs[0].set_title('time: %d a.u.  $E=%.4f$a.u.' % (t[i], energy[0,i]))
            axs[1].set_xlabel(r'$t$ / a.u.')

            plt.tight_layout()
            # for gif
            if gif:
                gif_frame(i)

            # end gif
            if plot_online:
                plt.show()
                plt.pause(pause)

########## reading input.q ##########
if exists('input.q'):
    namelist = f90nml.read('input.q')
    print('\nReading namelist in file input.q\n\n', namelist)
    run = namelist['general']['run']
    rank = namelist['general']['rank']
    ngrid = namelist['general']['ngrid']
    nstates = namelist['general']['nstates']
    if 'use_field' in namelist['general']:
        use_field = namelist['general']['use_field']
    else:
        use_field = False
else:
    print('File input.q not found. Following data required.')
    run = int(input('run: '))
    rank = int(input('rank: '))
    ngrid = int(input('ngrid: '))
    nstates = int(input('nstates: '))
    use_field = bool(input('use_field [True/False]: '))

########## initialize ##########
if run == 0:
    print("\nPlotting real time dynamics from Qdyn.\n")
elif run == 1:
    print("\nPlotting imaginary time dynamics from Qdyn.\n")

if rank == 3:
    print("ERROR: 3D plotting is not available! Exiting...")
    exit(1)

wf, en_states = [], []

wf, nframes_wf = read_wf()
en_states, nframes = read_energies()

if use_field:
    field = read_field()

# TODO: this should be changed when I add more states to RT propagation.
# It should be united with IM reading, I will just additionally read also energy.dat
if run == 0:
    energy = en_states[0][1:]
    t = en_states[0][0]

# estimate number of frames
if nframes != nframes_wf:
    print("Number of data in wf files and energy files is not matching. Exiting..")
    exit(1)

if final_plot and plot_online:
    plot_online = False

if not final_plot and not plot_online and not gif:
    print("*Either final_plot or plot_online or gif must be true.")
    exit(1)

if plot_online:
    print("*plotting on-the-fly")
    plt.ion()

if gif:
    fig_folder = 'frames'
    if not exists(fig_folder):
        makedirs(fig_folder)
    frames = []

########## code ##########

if rank == 1:
    # plot wf
    if run == 0:
        plot_rt_1d()
    elif run == 1:
        plot_it_1d()
elif rank == 2:
    if run == 0:
        if heatmap:
            plot_rt_2d_heatmap()
    elif run == 1:
        if heatmap:
            plot_it_2d_heatmap()
        else:
            plot_it_2d()

# if plot_on_the_fly:
if plot_online:
    plt.pause(10.0)
    plt.ioff()

if gif:
    imageio.mimsave('wavepacket.gif', frames, duration=duration)
    rmdir(fig_folder)
