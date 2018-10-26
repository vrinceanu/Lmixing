#
#  code to generate figure 2 in the manuscript
#
from Lmixing.rates import rate
import numpy as np

T = 10
rho = 100

nrange = np.arange(50,800,50)

with open("fig2.dat", "w") as file:
    file.write("# T = {0:.2g}, rho= {1:.2g}, L = 1\n#quantum\n".format(T,rho))
for n in nrange:
    k1 = rate(n,1,T,rho)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k1)
        )

with open("fig2.dat", "a") as file:
    file.write("\n#semiclassical\n")
nrange = np.arange(20,450,10)
for n in nrange:
    k2 = rate(n,1,T,rho,method='semiclassical')
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2)
        )

with open("fig2.dat", "a") as file:
    file.write("\n# P_and_S\n")
nrange = np.arange(20,450,10)
for n in nrange:
    k3 = rate(n,1,T,rho,method='P_and_S')
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k3)
        )

with open("fig2.dat", "a") as file:
    file.write("\n# T = {0:.2g}, rho= {1:.2g}, L = n-2\n#quantum\n".format(T,rho))
for n in nrange:
    k1 = rate(n,n-2,T,rho)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k1)
        )

with open("fig2.dat", "a") as file:
    file.write("\n#semiclassical\n")
nrange = np.arange(20,450,10)
for n in nrange:
    k2 = rate(n,n-2,T,rho,method='semiclassical')
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2)
        )

with open("fig2.dat", "a") as file:
    file.write("\n# P_and_S\n")
nrange = np.arange(20,450,10)
for n in nrange:
    k3 = rate(n,n-2,T,rho,method='P_and_S')
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k3)
        )
