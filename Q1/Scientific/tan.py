#! /usr/bin/python

from math import pi, tan

x = pi/20
for i in range(0,11):
    print("tan(%1.9f) = %4e \n" %(x*i, tan(x*i)))
