import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# making life easier
def l2hnorm(x, h):
    return np.linalg.norm(x)*np.sqrt(h)

# Colelctors for convergence analysis
err = []
hs = []

#Choose here if you want to see the plots for h = 0.2
plot = False

#Choose here if you want to use the asymetrical domain
asym = True


# Creating a matrix for the symetrical case
def make_A(n):
    L = np.diag((n-1)*[2],0)
    L += np.diag((n-2)*[-1],1)
    L += np.diag((n-2)*[-1],-1)
    return L

# Creating a similar matrix for the asymetrical case
def make_A_asym(n):
    L = np.diag((n-1)*[2.0],0)
    L += np.diag((n-2)*[-1.0],1)
    L += np.diag((n-2)*[-1.0],-1)
    L[0,0] = 1
    L[0,1] = -(2/3)
    return L


for k in range(0,4):
    # k refinement
    h = 0.2/(2**k)
    hs.append(h)
    n = int(1/h)

    # creating the domain and a finer domain for the exactg solution
    omega = np.linspace(0,1, n+1)
    omega2 = np.linspace(0,1, 100)
    if asym:
        omega = np.delete(omega, 1)
    
    # Calling the proper matrix creation
    A = make_A(n) if not asym else make_A_asym(n-1)

    L = A/(h**2)

    # exact solutions
    fn1 = lambda x: -1/2*np.power(x,2)+3/2*x+1
    fn2 = lambda x: -1*np.exp(x)+np.e*x+2

    # inner points of the computational domain
    wh = np.linspace(0+h, 1-h, n-1)
    if asym:
        wh = wh[1:]

    # creating the boundary conditions
    f1 = np.ones(wh.shape[0])
    f1[0] += 1/h/h if not asym else 1/3/h/h
    f1[-1] += 2/h/h

    f2 = np.exp(wh)
    f2[0] += 1/h/h if not asym else 1/3/h/h
    f2[-1] += 2/h/h

    # Solving the system
    u1 = np.linalg.solve(L.copy(),f1)
    u2 = np.linalg.solve(L.copy(),f2)

    # inserting boundary conditions
    u1 = np.insert(u1, 0, 1)
    u1 = np.append(u1, 2)

    u2 = np.insert(u2, 0, 1)
    u2 = np.append(u2, 2)

    # Plotting of the solution
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

    #Convergence analysis of u2
    print(l2hnorm(u1-fn1(omega), h))
    err.append(l2hnorm(u2-fn2(omega),h))
    print(err[-1])


#FInding alphha and C
to_fit = lambda x,a,c: c*np.power(x,a)

fitted = curve_fit(to_fit, hs, err)[0]
print(fitted)

# Plotting Convergence analysis
plt.loglog(hs,err,'-x', label="Global Error")
plt.loglog(hs, to_fit(hs, *fitted), label="Fitted logarithmic function")
plt.legend()
plt.axis('equal')
plt.show()