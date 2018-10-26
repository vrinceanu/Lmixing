"""
Gauss - Lobatto - Kronrod quadrature
reference: NR, Gautchi
"""
import mpmath

class quadGLK:
    def __init__(self, tol):
        self.tol = tol
        self.alpha = mpmath.sqrt(mpmath.mpf(2)/mpmath.mpf(3))
        self.beta = 1/mpmath.sqrt(mpmath.mpf(5))
        self.counter = 0
                                       

    def integrate(self, f, a, b):
        m = (a+b)/2
        h = (b-a)/2
        x1 = mpmath.mpf('.94288241569547971905635175843185720232')
        x2 = mpmath.mpf('.64185334234578130578123554132903188354')
        x3 = mpmath.mpf('.23638319966214988028222377349205292599')
        A = mpmath.mpf('.015827191973480183087169986733305510591')
        B = mpmath.mpf('.094273840218850045531282505077108171960')
        C = mpmath.mpf('.15507198733658539625363597980210298680')
        D = mpmath.mpf('.18882157396018245442000533937297167125')
        E = mpmath.mpf('.19977340522685852679206802206648840246')
        F = mpmath.mpf('.22492646533333952701601768799639508076')
        G = mpmath.mpf('.24261107190140773379964095790325635233')
        x = [-1, -x1, -self.alpha, -x2, -self.beta, -x3, 0, x3, self.beta, x2, self.alpha, x1, 1]
        y = [f(m + j*h) for j in  x]
        self.counter = 13
        fa = y[0]
        fb = y[12]
        i2 = h/6*(y[0] + y[12] + 5*(y[4] + y[8]))
        i1 = (h/1470)*(77*(y[0]+y[12])+432*(y[2]+y[10])+625*(y[4]+y[8])+672*y[6])
        ig = A*(y[0]+y[12]) + B*(y[1]+y[11]) + C*(y[2]+y[10]) + D*(y[3]+y[9]) + \
             E*(y[4]+y[8])  + F*(y[5]+y[7])  + G*y[6]
        err1 = abs(i1 - ig)
        err2 = abs(i2 - ig)
        ig =   abs(ig)
        return(self.adaptlob(f,a,b,fa,fb,ig))

    def adaptlob(self,f,a,b,fa,fb,ig):
        m = (a+b)/2
        h = (b-a)/2
        mll = m - self.alpha*h
        ml  = m - self.beta*h
        mr  = m + self.beta*h
        mrr = m + self.alpha*h
        fmll = f(mll)
        fml  = f(ml)
        fm   = f(m)
        fmr  = f(mr)
        fmrr = f(mrr)
        self.counter += 5
        i2  = h/6*(fa + fb + 5*(fml + fmr))
        i1  = h/1470*(77*(fa+fb) + 432*(fmll + fmrr) + 625*(fml + fmr) + 672*fm)
        if abs(i1-i2) <= self.tol*ig :
            return i1
        else:
            return(self.adaptlob(f, a,   mll, fa,   fmll, ig) + \
                   self.adaptlob(f, mll, ml,  fmll, fml,  ig) + \
                   self.adaptlob(f, ml,  m,   fml,  fm,   ig) + \
                   self.adaptlob(f, m,   mr,  fm,   fmr,  ig) + \
                   self.adaptlob(f, mr,  mrr, fmr,  fmrr, ig) + \
                   self.adaptlob(f, mrr, b,   fmrr, fb,   ig))
                   
if  __name__=="__main__":
    import mpmath

    mpmath.mp.dps = 100

    def f(x):
        return(mpmath.exp(-x))

    J = quadGLK(1.e-10)
    print(J.integrate(f, 0, 1)/(1 - 1/mpmath.e) - 1)
    print(J.counter)

    def f(x):
        return(mpmath.exp(-x*x))

    J = quadGLK(1.e-20); res = J.integrate(f, -20, 20)**2/mpmath.pi - 1
    print("# of evaluations: {0:4g}, error: {1:12.4e}".format(J.counter, float(res)))

    
