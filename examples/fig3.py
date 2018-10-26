#
#  code to generate figure 3 in the manuscript
#
import math
from Lmixing.rates import rate
import numpy as np

n_range = [14, 20, 29, 42, 60, 87, 126, 182, 264, 382, 554, 803, 1165]

x_range = [pow(10,-1 + j*0.5) for j in range(15)]
T_range = [math.sqrt(100*x) for x in x_range]

X,Y = np.meshgrid(x_range,n_range)

L = 1

Z1 = [[rate(n, L , T, rho = 100, method='quantum') for T in T_range] for n in n_range]]
Z2 = [[rate(n, L , T, rho = 100, method='semiclassical') for T in T_range] for n in n_range]
Z3 = [[rate(n, L , T, rho = 100, method='P_and_S') for T in T_range] for n in n_range]

np.savez('fig3', X, Y, Z1, Z2, Z3)
