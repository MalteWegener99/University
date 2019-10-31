import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from functools import partial
import time


def make_L(Nx, Ny):
    Dx = sp.diags((Nx-1)*[1])
    Dx += sp.diags((Nx-2)*[-1],-1)
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


def discretize(x_d, y_d, h):
    nx = int(x_d/h)
    ny = int(y_d/h)
    return make_L(nx,ny)/h/h

def get_grid(x_d, y_d, h):
    grid = np.mgrid[h:y_d:h, h:x_d:h]
    return (grid[1,:,:], grid[0,:,:])


def source(xx, yy, t):
    v = 4
    a = -50
    f = np.sin(2*np.pi*v*t)*np.exp(a*(np.square(xx-2*np.ones(xx.shape))+np.square(yy-2*np.ones(yy.shape))))
    return np.reshape(f, xx.shape[0]*xx.shape[1])

# Normal step
def step(u0, um1, A, dt, c, f):
    cdt = (c*dt)**2
    tmp = (2*sp.identity(A.shape[0])+cdt*(A))
    return tmp@u0+dt*dt*f-um1

# Initial step
def step0(u0, um1, A, dt, c, f):
    cdt = 0.5*(c*dt)**2
    tmp = (cdt*(A))
    return tmp@u0+dt*dt*f

# Unsteady solver interface
def unsteady_solver(u0, A, c, dt, T, saves, xx, yy):
    
    t = 0
    um1 = u0
    un = u0
    us = []
    if 0 in save_points:
        us.append(u0)
    while t <= T:
        start = time.time()
        f = source(xx, yy, t)
        if t == 0:
            ut = step0(un, um1, A, dt, c, f)
        else:
            ut = step(un, um1, A, dt, c, f)
        um1 = un
        un = ut
        t += dt
        for s in saves:
            if abs(s-t) <= dt/2:
                us.append(un)
        print("\rAt time %1.5f s, Last step took %2.8f s, expected time left: %3.2f s"%(t, time.time()-start, (T-t)/dt*(round(time.time()-start, 2))), end='')
    print("")
    return us

def initial(xx, yy):
    return np.reshape(np.zeros(xx.shape), xx.shape[0]*xx.shape[1])

c = 1.333
x = 4
y = 12
h = 0.02
dt = 0.99*h/c/(np.sqrt(2))

grid = get_grid(x,y,h)
reshaper = lambda u: np.reshape(u, [grid[0].shape[0], grid[0].shape[1]])[::-1,:]
L = -1 * discretize(x,y,h)

# Where to save the solution
save_points = [1, 2, 3, 4, 5, 6, 7, 300]

uno = unsteady_solver(initial(*grid), L, c, dt, 300 , save_points, *grid)

plt.imshow(reshaper(uno[-1]))
plt.colorbar()
plt.show()

print(len(uno))
fig = plt.figure()

mx = np.max(np.array(uno))
mn = np.min(np.array(uno))

# Normalize with 0 being at the middle
if abs(mx)>abs(mn):
    mn = -1*abs(mx)
else:
    mx = abs(mn)
# Plotting
for i in range(len(uno)):
    plt.subplot(180+i+1)
    plt.title("t = %1.1f"%(save_points[i]))
    im = plt.imshow(reshaper(uno[i]), vmax=mx, vmin=mn)#, cmap="gnuplot")

cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.1])
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')

plt.show()
