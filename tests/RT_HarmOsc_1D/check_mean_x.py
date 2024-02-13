import numpy as np

# inputs
grid_size = 512
x0 = 2 # initial position of the wave funciton (zero initial momentum considered)
exact_omega = 0.08 # exact omega from the model
thresh = 1e-4 # threshold for the error

# reading the wave function and reshaping
wf = np.genfromtxt('wf1d.1.out')
nframes_wf = int(np.shape(wf)[0] / grid_size)
wf = np.reshape(wf, (nframes_wf, grid_size, 5)).transpose((0, 2, 1))

# reading time from energy file
t = np.genfromtxt('energies.dat').T[0]

# calculating mean x position
x_mean = []
for i in range(0, nframes_wf):
    x_mean.append(np.trapz(y=wf[i, 3] * wf[i, 0], x=wf[i, 0]))

# calculating exact result
exact = x0 * np.cos(exact_omega * t)

# calculating cumulative error
error = np.sum((exact - x_mean) ** 2)

# printing result: True/False
print(error < thresh)
