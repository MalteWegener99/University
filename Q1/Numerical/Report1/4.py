import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la

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

def source(xx,yy):
    return 20*np.sin(np.pi*yy)*np.sin(1.5*np.pi*xx+np.pi)+30*yy

def sourcevec(xx,yy):
    return np.reshape(source(xx,yy), (xx.shape[0]*xx.shape[1]))


x = 2
y = 1
h = 0.01
grid = get_grid(x,y,h)
L = discretize(x,y,h)
solution = la.spsolve(L,sourcevec(*grid))
plt.imshow(np.reshape(solution, [grid[0].shape[0], grid[0].shape[1]]))
plt.show()