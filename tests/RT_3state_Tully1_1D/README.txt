"""Extended 3-state Tully model 1 created by Tomas Jira. Reference results are based on his simulations with ACRON.
The purpose of the test is to check numerical diagonalization of the Hamiltonian when building expH1 for more than 2 states where we
 don't have analytic formula.
.

The model is defined by the following Hamiltonian:
def pot(x):
    V=np.zeros((3,3))
    V[0,0]=0.01*np.tanh(0.5*x)
    V[1,1]=0
    V[2,2]=-V[0,0]
    V[0,1]=0.001*np.exp(-x**2)
    V[1,0]=V[0,1]
    V[1,2]=0.001*np.exp(-x**2)
    V[2,1]=V[1,2]
    V[0,2]=0
    V[2,0]=V[0,2]
    return(V)
"""