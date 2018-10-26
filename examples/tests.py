#
# example of using Lmixing library
#
from Lmixing.rates import rate

T = 10
ne = 10
L = 1

print("  n     quantum       semiclassical  PS64     \n"\
      "  ============================================")

for n in range(10,110,10):
    print("{0:4d}  {1:12.4e}  {2:12.4e}  {3:12.4e}".format(n,\
            rate(n,L,T,ne),\
            rate(n,L,T,ne,method='semiclassical'),\
            rate(n,L,T,ne,method='P_and_S')))

