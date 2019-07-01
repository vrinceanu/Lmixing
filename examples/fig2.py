#
#  code to generate figure 2 in the manuscript
#
from Lmixing.rates import rate
import numpy as np

T = 10
rho = 100

with open("fig2.dat", "a") as file:
    file.write("# T = {0:g} rho = {1:g}, L = 1\n".format(T,rho))

meth = "quantum"
nrange = range(50,50,50)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,1,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

meth = "semiclassical"
nrange = range(20,800,10)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,1,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

meth = "P_and_S"
nrange = range(20,800,10)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,1,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

meth = "PS-M"
nrange = range(20,800,10)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,1,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

with open("fig2.dat", "a") as file:
    file.write("# T = {0:g} rho = {1:g}, L = n-2 \n".format(T,rho))
    
meth = "quantum"
nrange = range(50,50,50)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,n-2,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

meth = "semiclassical"
nrange = range(20,800,10)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,n-2,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

meth = "P_and_S"
nrange = range(20,800,10)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,n-2,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))

meth = "PS-M"
nrange = range(20,800,10)
with open("fig2.dat", "a") as file:
    file.write("\n# "+meth+"\n")
for n in nrange:
    k2 = rate(n,n-2,T,rho,method=meth)
    with open("fig2.dat", "a") as file:
        file.write("{0:4d} {1:12.4e}\n".format(n,k2))