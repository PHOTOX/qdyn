"""Tully model 1 testing and comparing to results in ref https://doi.org/10.1021/acsomega.2c04843

potential defined as
AA=0.01
BB=0.6
CC=0.001
DD=1.0
E0=0.05

def pot(x):
    V=np.zeros((2,2))
    V[0,0]=AA*np.tanh(BB*x)
    V[1,1]=-V[0,0]
    V[0,1]=CC*np.exp(-DD*x*x)
    V[1,0]=V[0,1]
    return(V)
"""