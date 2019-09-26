import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from functools import partial
from numba import jit

def source(x,y):
    a = -5
    s = np.zeros(x.shape)
    for j in range(1,4,2):
        for i in range(1,8,2):
            print(i,j)
            s += np.exp(a*np.square(x-i*np.ones(x.shape))+a*np.square(y-j*np.ones(y.shape)))
    return s

def get_grid(x_d, y_d, h):
    grid = np.mgrid[h:y_d:h, h:x_d:h]
    return (grid[1,:,:], grid[0,:,:])

def k(x,y):
    return np.ones(x.shape)+4*x+6*y


def discretize(x_d, y_d, h, xx, yy):
    print("Start")
    nx = int(x_d/h)-1
    ny = int(y_d/h)-1
    h2x = np.ones(xx.shape)*h/2
    h2y = np.ones(yy.shape)*h/2
    kxp = np.reshape(k(xx+h2x, yy), nx*ny)
    kxm = np.reshape(k(xx-h2x, yy), nx*ny)
    kyp = np.reshape(k(xx, yy+h2y), nx*ny)
    kym = np.reshape(k(xx, yy-h2y), nx*ny)
    main_diag = (kxm.copy()+kxp.copy()+kyp.copy()+kym.copy())/h/h

    off_diag_up = -1*kxp.copy()[:-1:]/h/h
    t = off_diag_up[nx::nx]
    t[:] = 0

    off_diag_down = -1*kxm.copy()[1::]/h/h
    t = off_diag_down[nx::nx]
    t[:] = 0

    off_off_diag_up = -1*kyp.copy()[:(nx*ny-nx):]/h/h
    off_off_diag_down = -1*kym.copy()[nx::]/h/h

    L = sp.diags(main_diag)
    L += sp.diags(off_diag_up, 1)
    L += sp.diags(off_diag_up, -1)
    L += sp.diags(off_off_diag_up, nx)
    L += sp.diags(off_off_diag_down, -1*nx)
    return L

        




x = 8
y = 4
h = 0.01

grid = get_grid(x,y,h)
L = discretize(x,y,h, *grid)
src_v = np.reshape(source(*grid), [L.shape[0],1])
src = source(*grid)[::-1,:]
k_d = k(*grid)[::-1,:]
u = la.spsolve(L, src_v)
u = np.reshape(u, [grid[0].shape[0], grid[0].shape[1]])[::-1,:]
plt.imshow(u)
plt.colorbar()
plt.show()