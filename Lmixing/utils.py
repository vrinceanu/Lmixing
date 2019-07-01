import sys
import math
import mpmath
import fractions

def Gamma1(x):
    G = mpmath.quad(lambda t: mpmath.exp(-t)/t, [x, mpmath.inf])
    return(float(G))


def Gamma(x):
    return(float(mpmath.gammainc(0, x)))

def Harmonic(n):
    H = 0
    for j in range(1,n+1):
        H += 1/j
    return(H)

def Erf(x):
    return(float(mpmath.erf(x)))

MAX_FACT = 5000
# build a table of factorials ahead of time
F = [math.factorial(x) for x in range(0, MAX_FACT)]

# easy factorial with check
def fact(n):
    if (int(n) != n or n >= MAX_FACT or n < 0):
        print("*** triangular error ! factorial of ", n, file=sys.stderr)
        return 0
    else:
        return F[int(n)]

# the square of the big_delta factor
def Delta(a,b,c):
    return fractions.Fraction(fact(a+b-c)*fact(a-b+c)*fact(-a+b+c),
                              fact(a+b+c+1))

# the square of the wigner 6-j symbol
def Wigner6j(a,b,c,d,e,f):
    d1 = Delta(a,b,c)
    if (d1 == 0): return 0
    d2 = Delta(a,e,f)
    if (d2 == 0): return 0
    d3 = Delta(d,b,f)
    if (d2 == 0): return 0
    d4 = Delta(d,e,c)
    if (d2 == 0): return 0
    sum = 0
    k_min = int(max(a+b+c, a+e+f, d+b+f, d+e+c))
    k_max = int(min(a+b+d+e, a+c+d+f, b+c+e+f))
    for k in range(k_min, k_max + 1):
        sum += (-1)**k*fractions.Fraction(fact(k+1),
               fact(k-a-b-c)*fact(k-a-e-f)*fact(k-d-b-f)*fact(k-d-e-c)*
               fact(a+b+d+e-k)*fact(a+c+d+f-k)*fact(b+c+e+f-k))
    return d1*d2*d3*d4*sum*sum



class PQ:
    """
    A class for calculation of P(nl -> nl'; alpha)
    """
    def __init__(self, n, l1, l2):
        self.n = n
        self.W = []
        self.H = []
        self.counter = 0
        self.L_range = range(abs(l1-l2), 1 + min(l1+l2, n-1))
        if (l1 > n - 1 or l2 > n - 1):
            self.L_range = range(0)
        for L in self.L_range:
            self.W.append((2*L+1)*(2*l2+1)*
            fractions.Fraction(fact(n-L-1),fact(n+L))*
            pow(4,L)*Wigner6j(l2,l1,L,(n-1)/2,(n-1)/2,(n-1)/2))
            h = []
# calculates the coefficients of Gegenbauer polynomial
            for k in range(0, math.floor((n - L - 1)/2) + 1):
                h.append((-1)**k*2**(n - L - 1 - 2*k)*
                         fractions.Fraction(fact(n - 1 - k),
                         fact(k)*fact(n - L - 1 - 2*k)))
            self.H.append(h)

    def at(self, a):
        self.counter += 1
        mpmath.mp.dps = 600
        alpha = mpmath.mpf(a)
        z = (1 + alpha*alpha*mpmath.cos(mpmath.pi*mpmath.sqrt(1 + alpha*alpha)))/(1 + alpha*alpha)
        ZZ = [1]
        G = []
        for k in range(0, self.n):
            ZZ.append(z*ZZ[-1])
        for L in self.L_range:
            g = 0
            h = self.H[L-self.L_range[0]]
            for k in range(0, math.floor((self.n - L - 1)/2) + 1):
                g += h[k].numerator*ZZ[self.n-L-1-2*k]
            G.append(pow(1 - z*z, L)*g*g)
        pq = sum([x.numerator*y/x.denominator for (x,y) in zip(self.W, G)])
        return(pq)
