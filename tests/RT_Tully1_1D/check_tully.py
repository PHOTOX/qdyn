from os.path import exists

from numpy import genfromtxt, linspace, interp, std, abs

### inputs
pop_thresh = 1e-1 # threshold for the populations std
ethresh = 1e-5 # threshold for the energy std and drift

### checking populations
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
pop_std = std(pop_ad_int - pop_ref_int)

### checking energies
# reading energies
energy = genfromtxt('energies.dat').T

# calculating energy standard deviation
energy_std = std(energy[1])

# calculating energy drift in percent of initial total energy
energy_drift = abs((energy[1,0] - energy[1,-1])/energy[1,0])

### deciding result
if (pop_std < pop_thresh) and (energy_drift < ethresh) and (energy_std < ethresh):
    result=True
else:
    result=False

# printing result: True/False
print(result)
