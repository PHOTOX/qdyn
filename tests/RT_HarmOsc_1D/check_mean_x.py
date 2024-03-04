import numpy as np

# inputs
grid_size = 512
x0 = 2 # initial position of the wave funciton (zero initial momentum considered)
exact_omega = 0.08 # exact omega from the model
thresh = 1e-4 # threshold for the error
ethresh = 2e-5 # threshold for the energy std and drift 

# reading the wave function and reshaping
wf = np.genfromtxt('wf1d.1.out')
nframes_wf = int(np.shape(wf)[0] / grid_size)
wf = np.reshape(wf, (nframes_wf, grid_size, 5)).transpose((0, 2, 1))

# reading energy file
energy = np.genfromtxt('energies.dat').T

# setting time
t = energy[0]

# calculating mean x position
x_mean = []
for i in range(0, nframes_wf):
    x_mean.append(np.trapz(y=wf[i, 3] * wf[i, 0], x=wf[i, 0]))

# calculating exact result
exact = x0 * np.cos(exact_omega * t)

# calculating cumulative error
error = np.sum((exact - x_mean) ** 2)

# calculating energy standard deviation
energy_std = np.std(energy[1])

# calculating energy drift in percent of initial total energy
energy_drift = np.abs((energy[1,0] - energy[1,-1])/energy[1,0])

# deciding result
if (error < thresh) and (energy_drift < ethresh) and (energy_std < ethresh):
    result=True
else:
    result=False

# printing result: True/False
print(result)
