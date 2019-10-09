#! /usr/bin/python
from bigfloat import *

# Automatic big ints
x = 23**45
# Table of 23^45
for i in range(1,11):
    print("%d * %d = %d" %(x,i,x*i))

# couldnt get bigfloat to install so this is not tested
with precision(133):
    pi = BigFloat("3.14159265358979323846264338327950288419716939937510")
    pi2 = div(pi,2)
    for i in range(1,11):
        print("tan(%.40f) = %.40f"%(mul(pi2, i), tan(mul(pi2, 2))))
