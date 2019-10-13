import numpy as np
import matplotlib.pyplot as plt

poly = lambda x: x*x*x - 12*x*x + 23*x -10
der1 = lambda x: 3*x*x - 24*x + 23
der2 = lambda x: 6*x - 24

x = np.arange(-10,10,1)
y = poly(x)
yder = der1(x)
yder2 = der2(x)
yderfd = np.gradient(y, x)

fig = plt.figure()
plt.subplot(221)
plt.title("f(x)")
plt.plot(x,y)

plt.subplot(222)
plt.title("f'(x)")
plt.plot(x,yder)

plt.subplot(223)
plt.title("f''(x)")
plt.plot(x,yder2)

plt.subplot(224)
plt.title("approximately f(x)")
plt.plot(x,yderfd)

plt.show()