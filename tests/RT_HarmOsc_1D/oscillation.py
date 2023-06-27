import numpy as np
from scipy.optimize import curve_fit

def fit_cos(x, omega, phi, a, b): # Morse
    return a*np.cos(x*omega + phi)+b

e = np.genfromtxt('energies.dat').T

fit,cov = curve_fit(fit_cos, xdata=e[0], ydata=e[2], p0=(0.1,10,0.01,-1))

print(fit[0])