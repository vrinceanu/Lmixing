#
#  code to check the Gauss Lobatto Kronrod (GLK) quadrature
#
import mpmath
from Lmixing.utils import *
from Lmixing.quadGLK import *

n=10
L=1
prob = PQ(n,L,L+1)
Q = quadGLK(1.e-5)

theta = 1.e-2
mpmath.mp.dsp = 40

def W(z):
    return(prob.at(3/2/n/z)*z*mpmath.exp(-theta*z**2/2))

zmin = mpmath.mpf(1.e-10)
zmax = 10/mpmath.sqrt(theta)
res = Q.integrate(W, zmin, zmax)

m_result = 3.518265799
print("10,1 -> 10,2 integral for I(theta = 0.01)")
print("Lmixing result: {0:10.8f}   Mathematica result: {1:10.8f}\n".format(float(res),m_result))

n=500
L=498
prob = PQ(n,L,L+1)
Q = quadGLK(1.e-4)

theta = 1.e-2
mpmath.mp.dsp = 400

def W(z):
    return(prob.at(3/2/n/z)*z*mpmath.exp(-theta*z**2/2))

zmin = mpmath.mpf(1.e-10)
zmax = 10/mpmath.sqrt(theta)
res = Q.integrate(W, zmin, zmax)

m_result = 0.02741316138011
print("500,498 -> 500,499 integral for I(theta = 0.01)")
print("Lmixing result: {0:10.8f}   Mathematica result: {1:10.8f}\n".format(float(res),m_result))