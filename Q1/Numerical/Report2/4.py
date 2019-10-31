import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from functools import partial
import time
from mpl_toolkits.mplot3d import Axes3D  

# logging stuff
iterations = []
last_iter = {}
norms = []

def make_L(Nx, Ny):
    Dx = sp.diags((Nx-1)*[1.])
    Dx += sp.diags((Nx-2)*[-1.],-1)
    rowx = sp.csr_matrix((1,Nx-1))
    rowx[0,-1] = -1
    Dx = sp.vstack((Dx, rowx))
    Lx = Dx.transpose().dot(Dx)

    Dy = sp.diags((Ny-1)*[1])
    Dy += sp.diags((Ny-2)*[-1],-1)
    rowy = sp.csr_matrix((1,Ny-1))
    rowy[0,-1] = -1
    Dy = sp.vstack((Dy, rowy))
    Ly = Dy.transpose().dot(Dy)
    return sp.kronsum(Lx,Ly)

def rect_mask(xx, yy ,x1, x2, y1, y2):
    maskx1 = np.zeros(xx.shape, dtype=bool)
    maskx1[xx>=x1] = True
    maskx2 = np.zeros(xx.shape, dtype=bool)
    maskx2[xx<=x2] = True
    masky1 = np.zeros(xx.shape, dtype=bool)
    masky1[yy>=y1] = True
    masky2 = np.zeros(xx.shape, dtype=bool)
    masky2[yy<=y2] = True
    return np.logical_and(np.logical_and(maskx1, masky1), np.logical_and(maskx2, masky2))

    

def discretize(x_d, y_d, h):
    nx = int(x_d/h)
    ny = int(y_d/h)
    return make_L(nx,ny)/h/h

def get_grid(x_d, y_d, h):
    grid = np.mgrid[h:y_d:h, h:x_d:h]
    return (grid[1,:,:], grid[0,:,:])

def k(xx, yy):
    alpha = 10

    masker = lambda x1, x2, y1, y2: rect_mask(xx, yy, x1, x2, y1, y2)
    r1 = masker(1,2,1,2)
    r2 = masker(1,3,3,5)
    r3 = masker(4,7,4,7)
    r4 = masker(9,12,4,6)
    r5 = masker(13,15,1,3)

    # Combinate
    R = np.logical_or(r1,r2)
    R = np.logical_or(R,r3)
    R = np.logical_or(R,r4)
    R = np.logical_or(R,r5)

    return np.reshape(np.where(R, alpha*np.ones(xx.shape), np.zeros(xx.shape)), xx.shape[0]*xx.shape[1]) 


def newton_raphson(u0, A, k, dt, epsilon):
    ui = u0.copy()
    i = 0
    while True:
        i += 1
        jacobian = A + sp.eye(A.shape[0]) - 2*sp.diags(ui)
        v = la.spsolve(sp.eye(A.shape[0])-dt*jacobian, ui - u0 - dt*(A@ui+k*ui*(1-ui)))
        ui = ui - v
        norms.append(np.linalg.norm(np.abs(v)))
        iterations.append(i)
        if np.linalg.norm(np.abs(v))<epsilon:
            if i in last_iter:
                last_iter[i] += 1
            else:
                last_iter[i] = 1
            return ui

def picard(u0, A, k, dt, epsilon):
    f = lambda u: A@u+k*u*(1-u)
    ui = u0.copy()
    i = 0
    while np.linalg.norm(ui-dt*f(ui)-u0) > epsilon:
        i+=1
        ui = dt*f(ui)+u0
        norms.append(np.linalg.norm(ui-dt*f(ui)-u0))
        iterations.append(i)
    if i in last_iter:
        last_iter[i] += 1
    else:
        last_iter[i] = 1
    return ui


def step_FE(u0, A, k, dt):
    return u0+dt*A@u0+dt*k*u0*(1-u0)

def step_BE(u0, A, k, dt):
    if not use_picard:
        return newton_raphson(u0, A, k, dt, 0.001)
    else:
        return picard(u0, A, k, dt, 0.00001)

def stable_dt(h):
    return (h**2)/4

def unsteady_solver(u0, A, dt, T, saves, k, method="FE"):
    step = step_FE
    if method == "BE":
        step = step_BE
        if not use_picard:
            dt = 0.4
        else:
            dt = dt/2

    t = 0
    un = u0
    us = []
    if 0 in save_points:
        us.append(u0)
    global_start = time.time();
    while t < T:
        start = time.time()
        un = step(un, A, k, dt)
        t += dt
        for s in saves:
            if abs(s-t) <= dt/2:
                us.append(un)
        print("\rAt time %2.5f s, Last step took %2.8f s, expected time left: %3.0f s, total time: %3.2fs"%(t, time.time()-start, (time.time()-global_start)/t*T-(time.time()-global_start), time.time()-global_start), end='')
    print("")
    return us

def initial(xx, yy):
    f = np.exp(-2*(np.square(xx-1.5*np.ones(xx.shape))+np.square(yy-1.5*np.ones(yy.shape))))
    return np.reshape(f, xx.shape[0]*xx.shape[1])

#specify if to use Picard with BE, default is NR
use_picard = False
x = 16
y = 8
h = 0.04
stable_dt = (h**2)/4*0.99

grid = get_grid(x,y,h)
reshaper = lambda u: np.reshape(u, [grid[0].shape[0], grid[0].shape[1]])[::-1,:]
A = -1 * discretize(x,y,h)

save_points = [0, 1, 2, 3, 5, 10, 20, 30, 40]
#save_points = [0, 1, 1.33, 1.66, 2, 2.5, 3, 8, 40]

# Set method to either "FE" or "BE"
uno = unsteady_solver(initial(*grid), A, stable_dt, 40 , save_points, k(*grid), method="FE")

# Solution Plotting
fig = plt.figure()
mx = np.max(np.array(uno))
mn = np.min(np.array(uno))
for i in range(len(uno)):
    plt.subplot(3,3, 0+(i+1))
    plt.title("t = %2.1fs"%(save_points[i]))
    im = plt.imshow(reshaper(uno[i]))#, vmax=mx, vmin=mn)#, cmap="gnuplot")

cbar_ax = fig.add_axes([0.92, 0.05, 0.05, 0.9])
fig.colorbar(im, cax=cbar_ax)

plt.show()

#Iteration plotting
fig = plt.figure()
if len(norms) > 0:
    for i in range(max(last_iter.keys())):
        if i not in last_iter:
            last_iter[i] = 0

    plt.subplot(1,2,1)
    plt.scatter(iterations, norms)
    plt.ylim([0,1.1*max(norms)])
    plt.xlabel("Iteration number")
    plt.ylabel("Residual")
    plt.subplot(1,2,2)
    plt.xlabel("# of iterations for convergence")
    plt.ylabel("Occurences")
    plt.bar(*zip(*last_iter.items()))
    plt.show()