from os.path import exists

from numpy import genfromtxt, linspace, interp, std

# reading calculated data
input_file = 'pop_ad.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
pop_ad = genfromtxt(input_file).T

# reading reference from https://doi.org/10.1021/acsomega.2c04843 (I used their code to produce it)
input_file = 'reference_ad_pop.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
ref_pop = genfromtxt(input_file).T

# interpolating data to have the same grid
t = linspace(max([pop_ad[0, 0], ref_pop[0, 0]]), min([pop_ad[0, -1], ref_pop[0, -1]]), 1000)
pop_ad_int = interp(x=t, xp=pop_ad[0], fp=pop_ad[1])
pop_ref_int = interp(x=t, xp=ref_pop[0], fp=ref_pop[1])

# calculating standard deviation
std = std(pop_ad_int - pop_ref_int)

if std > 0.01:
    print(False)
else:
    print(True)
