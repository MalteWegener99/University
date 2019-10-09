import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from functools import partial

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
    return 20*np.sin(np.pi*yy)*np.sin(1.5*np.pi*xx+np.pi)

def sourcevec(xx,yy):
    return np.reshape(source(xx,yy), (xx.shape[0]*xx.shape[1]))

def make_boundary_vec(x, y, xx, yy, f):
    Bx0 = f[0](xx[0,:])
    Bx1 = f[2](xx[0,:])
    By0 = f[1](yy[:,0])
    By1 = f[3](yy[:,0])

    return (Bx0, Bx1, By0, By1)

def unitary(n, first):
    x = np.zeros([n,1])
    if first:
        x[0] = 1
    else:
        x[-1] = 1
    return x

def apply_boundary(x, y, xx, yy, h):
    f1 = lambda x: np.sin(0.5*np.pi*x)
    f2 = lambda y: np.sin(2*np.pi*y)
    f3 = lambda x: 0*x
    f4 = lambda y: np.sin(2*np.pi*y)
    bds = make_boundary_vec(x, y, xx, yy, [f1, f2, f3, f4])
    ny = xx.shape[0]
    nx = xx.shape[1]
    B = sp.kron(unitary(ny,True),np.reshape(bds[0],[nx,1]))
    B += sp.kron(unitary(ny,False),np.reshape(bds[1],[nx,1]))
    B += sp.kron(np.reshape(bds[2],[ny,1]), unitary(nx,True))
    B += sp.kron(np.reshape(bds[3],[ny,1]), unitary(nx,False))
    return B/h/h

imshow = partial(plt.imshow, interpolation='nearest')

x = 2
y = 1
h = 0.02
grid = get_grid(x,y,h)
L = discretize(x,y,h)
B = apply_boundary(2,1,*grid, h)
plt.spy(L, precision=0.001, marker='o', markersize=9*h)
plt.show()
sv = np.reshape(sourcevec(*grid), B.shape)
solution = la.spsolve(L,sv+B)
imshow(source(*grid)[::-1,:])
plt.show()
imshow(np.reshape(solution, [grid[0].shape[0], grid[0].shape[1]])[::-1,:])
plt.colorbar()
plt.show()