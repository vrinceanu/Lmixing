import sys
import math
import mpmath
from .constants import *
from .utils import PQ, Gamma, Gamma1, Erf, Harmonic
from .quadGLK import quadGLK

'''
dimensional factors and maximum order for power expansion
'''
A_factor = (m_e/m_p)*(k_B/hartree)**2*(0.01/a0)**3/(4*pi)
unit_factor = math.sqrt(8*math.pi)*a0**2*v0*math.sqrt(m_p*v0**2/(k_B))/cm**3
MAX_ORDER = 16

def kQ(n,li,lf,theta):
    prob = PQ(n, li, lf)
    W = lambda z: prob.at(3/2/n/z)*z*mpmath.exp(-theta*z*z/2)
    tol = 1.e-4
    Q = quadGLK(tol)
    zmax = max(math.sqrt(30./theta), 20)
    k1  = Q.integrate(W, 1.e-8, 0.1)
    k1 += Q.integrate(W, 0.1, 2)
    k1 += Q.integrate(W, 2, 10)
    k1 += Q.integrate(W, 10, zmax)
    return(float(k1))

def kQd(n,l,theta):
    prob1 = PQ(n, l, l+1)
    prob2 = PQ(n, l, l-1)
    tol = 1.e-4
    zmax = max(math.sqrt(30./theta), 20)
    W = lambda z: (prob1.at(3/2/n/z) + prob2.at(3/2/n/z))*z*mpmath.exp(-theta*z*z/2)
    Q = quadGLK(tol)
    k1  = Q.integrate(W, 1.e-8, 0.1)
    k1 += Q.integrate(W, 0.1, 2)
    k1 += Q.integrate(W, 2, 10)
    k1 += Q.integrate(W, 10, zmax)
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

def kPS_M(n,l,theta):
    Dnl = 6*(1 - l**2/n**2 - l/n**2 - 1/n**2)
    P1 = 1/4;
    beta = Dnl*theta/4/P1
    k = ((math.sqrt(math.pi)/2*Erf(math.sqrt(beta))/math.sqrt(beta) - math.exp(-beta))/beta + Gamma(beta))*Dnl/4
    return(k)

def kPS(n,l,theta):
    Dnl = 6*(1 - l**2/n**2 - l/n**2 - 1/n**2)
    k = Dnl/4*(1 - EulerGamma - math.log(Dnl*theta/2))
    return(k)


def rate(n, li, T, rho = 0, lf = None, method = 'quantum', mu = 1/2):
    """
    function to calculate L-mixing rate coefficients in cm^3/s with various approximations
    required arguments:
        n   -- principal quantum number
        li  -- angular momentum quantum number for initial state
        T   -- proton temperature in K
    keyword arguments:
        rho -- electron density in cm^-3 (default: 0)
        lf  -- final angular momentum quantum number (default: not used)
        method -- "quantum" (default), "semiclassical", "classical", "P_and_S", "Born", and "PS-M"
        mu  -- reduced mass / proton mass (default is mu = 0.5 for H+ collisions with H(n,l))
    when lf is not provided, rho is required and the rate is for the combined n li -> n li +/- 1 transitions
    when lf is given such that abs(lf - li) > 1, rho is not required and only the classical approximation is available
    """
    k = 0

    choices = ("quantum", "classical", "semiclassical", "Born", "P_and_S", "PS-M")
    if method not in choices:
        print("Posible methods are: ", choices)
        return(k)


    if (lf == None or abs(lf-li) == 1) and rho == 0:
        print("*** rho can not be zero for a dipole allowed transition ***",
              file=sys.stderr)
        return(k)

    theta = rho*n**4/T**2/A_factor*mu
    factor = unit_factor*n**4/math.sqrt(T)*math.sqrt(mu)
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
    if method=="PS-M":
        if lf:
            print(('*** Lf = {0}, {1} approximation works only '
                  'for dipole allowed transitions ***').format(lf,method),file=sys.stderr)
        else:
            k = kPS_M(n,li, theta)

    return(factor*k)
