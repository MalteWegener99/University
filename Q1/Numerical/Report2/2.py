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


def FE_step(u0, A, dt):
    return (sp.identity(A.shape[0])+(A)*dt)._mul_vector(u0)

def BE_step(u0, A, dt):
    return la.spsolve((sp.identity(A.shape[0])-(A)*dt), u0)

def stable_dt(h):
    return (h**2)/4

# Unsterady solver interface
def unsteady_solver(u0, A, dt, T, saves, method="FE", stabilizer=None):
    stepper = FE_step
    if method == "BE":
        stepper = BE_step
    if stabilizer is not None and method == "FE":
        dt = stabilizer()
    
    t = 0
    un = u0
    us = []
    if 0 in save_points:
        us.append(u0)
    while t <= T:
        t += dt
        un = stepper(un, A, dt)
        for s in saves:
            if abs(s-t) <= dt/2:
                us.append(un)

    return us

def initial(xx, yy):
    a = -5
    return np.reshape(np.exp(a*np.ones(xx.shape)*(np.square(xx-2*np.ones(xx.shape))+np.square(yy-2*np.ones(yy.shape)))), xx.shape[0]*xx.shape[1])

# Domain is [0,x] x [0,4]
x = 4
y = 4
h = 0.02
# specify  dt
dt = 0.015

grid = get_grid(x,y,h)
reshaper = lambda u: np.reshape(u, [grid[0].shape[0], grid[0].shape[1]])[::-1,:]
L = -1 * discretize(x,y,h)

# Specify where to save the solution
save_points = [0, 0.045, 0.09, 0.15]

print(stable_dt(h))

start = time.time()
# Remove the stabilizer arguemnt if you want to use dt as the timestep for FE
# Set method to either "FE" or "BE"
uno = unsteady_solver(initial(*grid), L, dt, 0.15, save_points, method="FE", stabilizer = lambda : stable_dt(h))
print(time.time()-start)
start = time.time()
dos = unsteady_solver(initial(*grid), L, dt, 0.15, save_points, method="BE")
print(time.time()-start)

# Normalization
mx = max(np.max(np.array(uno)), np.max(np.array(dos)))
mn = max(np.min(np.array(uno)), np.min(np.array(dos)))

# Plotting
fig = plt.figure()
plt.subplot(2,5,1)
plt.imshow(reshaper(uno[0]), vmin=mn, vmax=mx)
plt.title("t = 0s")
plt.subplot(2,5,2)
plt.imshow(reshaper(uno[1]), vmin=mn, vmax=mx)
plt.title("t = 0.045s")
plt.subplot(2,5,3)
plt.imshow(reshaper(uno[2]), vmin=mn, vmax=mx)
plt.title("t = 0.09s")
plt.subplot(2,5,4)
plt.imshow(reshaper(uno[3]), vmin=mn,  vmax=mx)
plt.title("t = 0.15s")

plt.subplot(2,5,6)
plt.imshow(reshaper(dos[0]), vmin=mn, vmax=mx)
plt.subplot(2,5,7)
plt.imshow(reshaper(dos[1]), vmin=mn, vmax=mx)
plt.subplot(2,5,8)
plt.imshow(reshaper(dos[2]), vmin=mn, vmax=mx)
plt.subplot(2,5,9)
im = plt.imshow(reshaper(dos[3]), vmin=mn, vmax=mx)

cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.show()
