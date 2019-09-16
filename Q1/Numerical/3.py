import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def l2hnorm(x):
    return np.sqrt(x.shape[0])*np.linalg.norm(x)

err = []
hs = []
plot = True

def make_A(n):
    L = np.diag((n-1)*[2],0)
    L += np.diag((n-2)*[-1],1)
    L += np.diag((n-2)*[-1],-1)
    return L

for k in range(0,5):
    h = 0.2/(2**k)
    hs.append(h)
    n = int(1/h)

    omega = np.linspace(0,1, n+1)

    A = make_A(n)

    L = A/(h**2)


    fn1 = lambda x: -1/2*np.power(x,2)+3/2*x+1
    fn2 = lambda x: -1*np.exp(x)+np.e*x+2

    wh = np.linspace(0+h, 1-h, n-1)
    f1 = np.ones(n-1)
    f1[0] += 1/h/h
    f1[-1] += 2/h/h

    f2 = np.exp(wh)
    f2[0] = 1/h/h
    f2[-1] += 2/h/h

    u1 = np.linalg.solve(L.copy(),f1)
    u2 = np.linalg.solve(L.copy(),f2)

    u1 = np.insert(u1, 0, 1)
    u1 = np.append(u1, 2)

    u2 = np.insert(u2, 0, 1)
    u2 = np.append(u2, 2)

    if k==0 and plot:
        plt.plot(omega, u1, label="numerical")
        plt.plot(omega, fn1(omega), label="analytical")
        plt.legend()
        plt.show()

        plt.plot(omega, u2, label="numerical")
        plt.plot(omega, fn2(omega), label="analytical")
        plt.legend()
        plt.show()

    # print(l2hnorm(u1-fn1(omega)))
    print(np.linalg.norm(u2-fn2(omega)))
    err.append(np.linalg.norm(u2-fn2(omega))*np.sqrt(n))

to_fit = lambda x,a,c: c*np.power(x,a)

plt.loglog(hs,err,'-x')
fitted = curve_fit(to_fit, hs, err)[0]
print(fitted)
plt.loglog(hs, to_fit(hs, *fitted))
plt.axis('equal')
plt.show()