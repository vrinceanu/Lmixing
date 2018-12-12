#
# example of using Lmixing library
#
from Lmixing.rates import rate
import time

'''
print("timing information")
start = time.time()
n = 200; L = 10; T = 10; ne = 100
r = rate(n,L,T,ne)
end = time.time()
print("rate coeff for n=",n," l=",L," T=",T,"K ne=",ne, " cm^-3")
print("is ", r," cm^3/2 and took ", end-start," sec")

start = time.time()
n = 200; L = 100; T = 10; ne = 100
r = rate(n,L,T,ne)
end = time.time()
print("rate coeff for n=",n," l=",L," T=",T,"K ne=",ne, " cm^-3")
print("is ", r," cm^3/2 and took ", end-start," sec")

start = time.time()
n = 400; L = 100; T = 10; ne = 100
r = rate(n,L,T,ne)
end = time.time()
print("rate coeff for n=",n," l=",L," T=",T,"K ne=",ne, " cm^-3")
print("is ", r," cm^3/2 and took ", end-start," sec")

start = time.time()
n = 400; L = 200; T = 10; ne = 100
r = rate(n,L,T,ne)
end = time.time()
print("rate coeff for n=",n," l=",L," T=",T,"K ne=",ne, " cm^-3")
print("is ", r," cm^3/2 and took ", end-start," sec")
'''

T = 10; ne = 10; L = 1

print("  n   quantum    semiclassical PS64     Born\n"\
      "  ================================================")

for n in range(10,50,10):
    print("{0:4d}  {1:10.4e} {2:10.4e}  {3:10.4e} {4:10.4e}".
    format(n,\
            rate(n, L, T, ne),\
            rate(n, L, T, ne, method='semiclassical'),\
            rate(n, L, T, ne, method='P_and_S'),\
            rate(n, L, T, ne, method='Born')
            ))
print( rate(20, 1, 10, rho=10))
            
print("\n  Delta L  quantum    classical\n"\
        "  ================================")
for h in range(2,8):
    print("{0:4d}       {1:10.4e} {2:10.4e}".format(h,
        rate(20, L, T, lf = L + h),
        rate(20, L, T, lf = L + h, method="classical")
      ))
