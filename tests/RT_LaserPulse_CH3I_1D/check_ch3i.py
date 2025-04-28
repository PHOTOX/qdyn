from os.path import exists

from numpy import genfromtxt, linspace, interp, std, trapz, sqrt, array, reshape, shape

xngrid = 256

# reading calculated data
wf = genfromtxt('wf1d_ad.1.dat')
nframes_wf = int(shape(wf)[0] / xngrid)
wf = reshape(wf, (nframes_wf, xngrid, 5)).transpose((0, 2, 1))
energy = genfromtxt('energies.dat').T
t = energy[0]

input_file = 'pop_ad.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
pop_ad = genfromtxt(input_file).T

# reading reference from Basile's & Fede's code for QD and exact factorization
input_file = 'reference_ad_pop.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
ref_pop = genfromtxt(input_file).T
input_file = 'ref_aver_x.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
ref_x = genfromtxt(input_file).T

# calculating mean <x> and <delta x>
x = [[], [], []]
for i in range(0, nframes_wf, 10):
    x[0].append(t[i] * 0.02418884254)
    dx = wf[i][0][1] - wf[i][0][1]
    norm = trapz(wf[i][3], wf[i][0], dx=dx)
    aver_x = trapz(wf[i][3] * wf[i][0], wf[i][0], dx=dx) / norm
    x[1].append(aver_x)
    aver_dx = sqrt(trapz(wf[i][3] * wf[i][0] ** 2, wf[i][0], dx=dx) / norm - aver_x ** 2)
    x[2].append(aver_dx)
x = array(x)

# interpolating data to have the same grid
t_int = linspace(max([pop_ad[0, 0], ref_pop[0, 0]]), min([pop_ad[0, -1], ref_pop[0, -1]]), 1000)
pop_ad_int = interp(x=t_int, xp=pop_ad[0], fp=pop_ad[1])
x_int = interp(x=t_int, xp=x[0], fp=x[1])
dx_int = interp(x=t_int, xp=x[0], fp=x[2])
pop_ref_int = interp(x=t_int, xp=ref_pop[0], fp=ref_pop[1])
x_ref_int = interp(x=t_int, xp=ref_x[0], fp=ref_x[1])
dx_ref_int = interp(x=t_int, xp=ref_x[0], fp=ref_x[2])

# calculating standard deviation
pop_std = std(pop_ad_int - pop_ref_int)
x_std = std(x_int - x_ref_int)
dx_std = std(dx_int - dx_ref_int)

if pop_std > 1e-3 or x_std > 5e-4 or dx_std > 5e-4:
    print(False)
else:
    print(True)
