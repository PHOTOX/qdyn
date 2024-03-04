"""Double well potential tested to reference code in ref https://doi.org/10.1021/acsomega.2c04843, which can be downloaded from zenodo

AA=0.001
BB=0.001
CC=0.001
DD=1.0

V = [[], []]
V_0_0 = AA*(x)**2
V_1_1 = AA*(x - 5)**2 + BB
V_0_1 = CC*np.exp(-DD*(x - 2.6)**2)
V_1_0 = V_0_1
"""