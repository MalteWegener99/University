import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def l2hnorm(x):
    return np.linalg.norm(x)/np.sqrt(x.shape[0])

err = []
hs = []
plot = True
asym = True

def make_A(n):
    L = np.diag((n-1)*[2],0)
    L += np.diag((n-2)*[-1],1)
    L += np.diag((n-2)*[-1],-1)
    return L

def make_A_asym(n):
    L = np.diag((n-1)*[2.0],0)
    L += np.diag((n-2)*[-1.0],1)
    L += np.diag((n-2)*[-1.0],-1)
    L[0,0] = 1
    L[0,1] = -(2/3)
    return L


for k in range(0,4):
    h = 0.2/(2**k)
    hs.append(h)
    n = int(1/h)

    omega = np.linspace(0,1, n+1)
    omega2 = np.linspace(0,1, 100)
    if asym:
        omega = np.delete(omega, 1)

    A = make_A(n) if not asym else make_A_asym(n-1)

    L = A/(h**2)


    fn1 = lambda x: -1/2*np.power(x,2)+3/2*x+1
    fn2 = lambda x: -1*np.exp(x)+np.e*x+2

    wh = np.linspace(0+h, 1-h, n-1)
    if asym:
        wh = wh[1:]
    f1 = np.ones(wh.shape[0])
    f1[0] += 1/h/h if not asym else 1/3/h/h
    f1[-1] += 2/h/h

    f2 = np.exp(wh)
    f2[0] += 1/h/h if not asym else 1/3/h/h
    f2[-1] += 2/h/h

    u1 = np.linalg.solve(L.copy(),f1)
    u2 = np.linalg.solve(L.copy(),f2)

    u1 = np.insert(u1, 0, 1)
    u1 = np.append(u1, 2)

    u2 = np.insert(u2, 0, 1)
    u2 = np.append(u2, 2)

    if k==0 and plot:

        print(L)
        print(f1)
        print(f2)
        plt.plot(omega, u1, '-x', label="numerical")
        plt.plot(omega2, fn1(omega2), label="analytical")
        plt.legend()
        plt.show()

        plt.plot(omega, u2, '-x', label="numerical")
        plt.plot(omega2, fn2(omega2), label="analytical")
        plt.legend()
        plt.show()
        print(u1)
        print(fn1(omega))
        print(u2)
        print(fn2(omega))

    print(l2hnorm(u1-fn1(omega)))
    err.append(l2hnorm(u2-fn2(omega)))
    print(err[-1])
    input()

to_fit = lambda x,a,c: c*np.power(x,a)

fitted = curve_fit(to_fit, hs, err)[0]
print(fitted)

plt.loglog(hs,err,'-x', label="Global Error")
plt.loglog(hs, to_fit(hs, *fitted), label="Fitted logarithmic function")
plt.legend()
plt.axis('equal')
plt.show()