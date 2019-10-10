import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from functools import partial
import time
from mpl_toolkits.mplot3d import Axes3D  


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
    alpha = 1

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
    ui = u0
    while True:
        jacobian = A + sp.eye(A.shape[0]) - 2*sp.diags(ui)
        v = la.spsolve(sp.eye(A.shape[0])-dt*jacobian, ui - u0 - dt*(A@ui+k*ui*(1-ui)))
        ui = ui - v
        if((np.abs(v)<epsilon).all()):
            return ui


def step_FE(u0, A, k, dt):
    return u0+dt*A@u0+dt*k*u0*(1-u0)

def step_BE(u0, A, k, dt):
    return newton_raphson(u0, A, k, dt, 0.001)

def stable_dt(h):
    return (h**2)/4

def unsteady_solver(u0, A, dt, T, saves, k, method="FE"):
    step = step_FE
    if method == "BE":
        step = step_BE
        dt *= 10
    t = 0
    un = u0
    us = []
    if 0 in save_points:
        us.append(u0)
    while t <= T:
        start = time.time()
        un = step(un, A, k, dt)
        t += dt
        for s in saves:
            if abs(s-t) <= dt/2:
                us.append(un)
        print("\rAt time %1.5f s, Last step took %2.8f s, expected time left: %3.2f s"%(t, time.time()-start, (T-t)/dt*(round(time.time()-start, 5))), end='')
    print("")
    return us

def initial(xx, yy):
    f = np.exp(-2*(np.square(xx-1.5*np.ones(xx.shape))+np.square(yy-1.5*np.ones(yy.shape))))
    return np.reshape(f, xx.shape[0]*xx.shape[1])

x = 16
y = 8
h = 0.1
dt = (h**2)/4*0.99

grid = get_grid(x,y,h)
reshaper = lambda u: np.reshape(u, [grid[0].shape[0], grid[0].shape[1]])[::-1,:]
A = -1 * discretize(x,y,h)

save_points = [0, 1, 2, 3, 5, 10, 20, 30, 40]

#stabilizer = lambda : stable_dt(h)

uno = unsteady_solver(initial(*grid), A, dt, 40 , save_points, k(*grid), method="FE")

print(len(uno))
plt.figure()
mx = np.max(np.array(uno))
mn = np.min(np.array(uno))
for i in range(len(uno)):
    print(i)
    plt.subplot(3,3, 0+(i+1))
    plt.imshow(reshaper(uno[i]))#, vmax=mx, vmin=mn)#, cmap="gnuplot")

plt.show()