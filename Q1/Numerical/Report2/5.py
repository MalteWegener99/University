import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from functools import partial
import time

# Creating the MAtrix as describve din the report
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

    #Kronsum is jsut a cleaner way than creating Identity matrices and stuff
    return sp.kronsum(Lx,Ly)


def discretize(x_d, y_d, h):
    nx = int(x_d/h)
    ny = int(y_d/h)
    return make_L(nx,ny)/h/h

def get_grid(x_d, y_d, h):
    grid = np.mgrid[h:y_d:h, h:x_d:h]
    return (grid[1,:,:], grid[0,:,:])

def source(xx,yy):
    a = -10
    b = 5
    return np.exp(a*(np.square(xx-b)+np.square(yy-b)))

def sourcevec(xx,yy):
    return np.reshape(source(xx,yy), (xx.shape[0]*xx.shape[1]))



imshow = partial(plt.imshow)


#Domain is (0,0)x(x,y)
x = 10
y = 10

#choose grid spacing
h = 0.1

gamma = [-40, 0, 40]
solutions = []
sources = []
residues = []

for i in range(3):
    # Creating grid, L 
    grid = get_grid(x,y,h)
    L = discretize(x,y,h)
    L = L - gamma[i] * sp.eye(L.shape[0])


    #Creation of the source vector
    sv = sourcevec(*grid)

    #Solving the system
    residuals = []

    def cb(rk):
        print("\rIn iteration number %4d, rk is %1.5e"%(len(residuals)+1,rk), end="")
        residuals.append(rk)
    start = time.time()
    solution, succ = la.gmres(L, sv, maxiter=5000, restart=5000, tol=1e-12, callback=cb)
    # check if GMRES COnverged
    if succ == 0:
        print("\nGMRES Converged")
    elif succ > 0:
        print("\nGMRES Converged but given tolerance not achieved or max iterations reached")
    else:
        print("\nYeah, you made an oopsie")

    print("GMRES took %3.2fs"%(time.time()-start))
    residues.append(residuals)
    solutions.append(solution)
    sources.append(sv)

    print("This should be small:", np.linalg.norm(sv-L@solution)/np.linalg.norm(sv)-residuals[-1])

    start = time.time()
    solution = la.spsolve(L, sv)
    print("spsolve took %3.2fs"%(time.time()-start))

#Showing source then Solution
reshaper = lambda u: np.reshape(u, [grid[0].shape[0], grid[0].shape[1]])[::-1,:]
for i in range(3):
    plt.semilogy(residues[i], label="$\gamma=$%d"%(gamma[i]))
plt.legend()
plt.ylabel("Residual")
plt.xlabel("Iteration")
plt.show()


fig = plt.figure()
plt.subplot(2,2,1)
plt.title("Source function")
plt.imshow(reshaper(sources[0]))
plt.colorbar()
for i in range(3):
    plt.subplot(2,2,i+2)
    plt.title("$\gamma=$%d"%(gamma[i]))
    plt.imshow(reshaper(solutions[i]))
    plt.colorbar()

plt.show()