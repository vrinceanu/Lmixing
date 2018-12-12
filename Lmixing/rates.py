import sys
import math
import mpmath
from .constants import *
from .utils import PQ, Gamma, Erf, Harmonic
from .quadGLK import quadGLK

mu = m_p/2
A_factor = (m_e/mu)*(k_B/hartree)**2*(0.01/a0)**3/(4*pi)
unit_factor = math.sqrt(8*math.pi)*a0**2*v0*math.sqrt(mu*v0**2/(k_B))/cm**3
MAX_ORDER = 16

def kQ(n,li,lf,theta):
    prob = PQ(n, li, lf)
    W = lambda z: prob.at(3/2/n/z)*z*mpmath.exp(-theta*z*z/2)
    tol = 1.e-4
    Q = quadGLK(tol)
    k1  = Q.integrate(W, 1.e-8, 0.1)
    k1 += Q.integrate(W, 0.1, 2)
    k1 += Q.integrate(W, 2, 100)
    return(float(k1))

def kQd(n,l,theta):
    prob1 = PQ(n, l, l+1)
    prob2 = PQ(n, l, l-1)
    zmin = mpmath.mpf(1.e-10)
    zmax = 6/math.sqrt(theta)
    tol = 1.e-4
    W = lambda z: (prob1.at(3/2/n/z) + prob2.at(3/2/n/z))*z*mpmath.exp(-theta*z*z/2)
    Q = quadGLK(tol)
    k1 = Q.integrate(W,zmin, zmax)
    return(float(k1))

def kC(n,li,lf):
    return(3/4*((li+lf) - min(li,lf)**2*(li+lf+2*abs(li-lf))/n**2)/(li+1/2)/abs(li-lf)**3)

def kSC(n,l,theta):
    GS = 0
    Dnl = 6*(1 - l**2/n**2 - l/n**2 - 1/n**2)
    x = 3*Dnl*theta/4
    eta = 0.2778545416190065
    for j in range(MAX_ORDER+1):
        A = 4**j/(j+2)/(j+3)/(2*j+3)/math.factorial(j)/math.factorial(2*j+1)
        B = 1/(j+2) + 1/(j+3) + 2/(2*j+3) + Harmonic(j) + 2*Harmonic(2*j+1)
        GS += A*(B - 3*EulerGamma - math.log(4*x))*x**j
    J = -math.exp(-eta**2*x)*3*eta/2/x
    J += 3*math.sqrt(pi)*Erf(eta*math.sqrt(x))/4/pow(x, 1.5)
    k = Dnl*(J + 18*GS)/4
    return(k)

def kB(n,l,theta):
    Dnl = 6*(1 - l**2/n**2 - l/n**2 - 1/n**2)
    k = (1 - math.exp(-Dnl*theta/2))/(2*theta) + Dnl*Gamma(Dnl*theta/2)/4
    return(k)

def kPS(n,l,theta):
    Dnl = 6*(1 - l**2/n**2 - l/n**2 - 1/n**2)
    k = Dnl/4*(1 - EulerGamma - math.log(Dnl*theta/2))
    return(k)


def rate(n, li, T, rho = 0, lf = None, method = 'quantum'):
    """
    function to calculate L-mixing rate coefficients in cm^3/s with various approximations
    required arguments:
        n   -- principal quantum number
        li  -- angular momentum quantum number for initial state
        T   -- proton temperature in K
    keyword arguments:
        rho -- electron density in cm^-3 (default: 0)
        lf  -- final angular momentum quantum number (default: not used)
        method -- "quantum" (default), "semiclassical", "classical", "P_and_S", "Born"
    when lf is not provided, rho is required and the rate is for the combined n li -> n li +/- 1 transition
    when lf is given such that abs(lf - li) > 1, rho is not required and only the classical approximation is available
    """
    k = 0
    if (lf == None or abs(lf-li) == 1) and rho == 0:
        print("*** rho can not be zero for a dipole allowed transition ***",
              file=sys.stderr)
        return(k)

    theta = rho*n**4/T**2/A_factor
    factor = unit_factor*n**4/math.sqrt(T)

    if method=="quantum":
        if lf:
            k = kQ(n,li,lf,theta)
        else:
            k = kQd(n,li,theta)
    if method=="classical":
        if lf and (abs(lf-li) > 1):
            k = kC(n,li,lf)
        else:
           print("*** Delta L has to be > 1 for classical approximation ***", file=sys.stderr)
    if method=="semiclassical":
        if lf:
            print(('*** Lf = {0}, {1} approximation works only '
                  'for dipole allowed transitions ***').format(lf,method),file=sys.stderr)
        else:
            k = kSC(n,li, theta)
    if method=="Born":
        if lf:
            print(('*** Lf = {0}, {1} approximation works only '
                  'for dipole allowed transitions ***').format(lf,method),file=sys.stderr)
        else:
            k = kB(n,li, theta)
    if method=="P_and_S":
        if lf:
            print(('*** Lf = {0}, {1} approximation works only '
                  'for dipole allowed transitions ***').format(lf,method),file=sys.stderr)
        else:
            k = kPS(n,li, theta)
    return(factor*k)
