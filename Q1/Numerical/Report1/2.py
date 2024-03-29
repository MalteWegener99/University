import numpy as np
import matplotlib.pyplot as plt

n = 10
domain = np.linspace(0,1,n+1)
h = 1/n
# Build L
L = np.diag((n-1)*[2],0)
L += np.diag((n-2)*[-1],1)
L += np.diag((n-2)*[-1],-1)
L = L/(h**2)

# Printign for Sanitycheck
print(L[0,:])
print(L[1,:])
print(L[-1,:])

# Show the matrix structure
plt.spy(L, marker="o", color='r')
plt.show()

#Calculate the eingenvalues
numerical_eigenvalues = [x for x in np.linalg.eigh(L)[0]]
analyticl_eigenvalues_dc = [(2/h)**2*np.sin(np.pi*i/2/n)**2 for i in range(1,n)]
analyticl_eigenvalues_op = [(np.pi * i)**2 for i in range(1,n)]

#show eigenvalues
plt.scatter([x.real for x in numerical_eigenvalues], [x.imag for x in numerical_eigenvalues], label="numerical discrete")
plt.scatter(analyticl_eigenvalues_dc, (n-1)*[0], label="analytical discrete", marker='x')
plt.scatter(analyticl_eigenvalues_op, (n-1)*[0], label="analytical operator")
plt.grid()
plt.legend()
plt.xlabel("Re")
plt.ylabel("Im")
print(analyticl_eigenvalues_dc)
print(analyticl_eigenvalues_op)
plt.show()


# Plotting the eigenvectors
for i in range(n-1):
    l = ([np.sin((i+1)*np.pi*j) for j in domain])
    plt.plot(domain,l, "+", markersize=20)


# Plotting the eigenfunctions
domain_f = np.linspace(0,1,10*n)
for i in range(1,n):
    print(i)
    plt.plot(domain_f, [np.sin(np.pi*i*x) for x in domain_f])
plt.show()


# Plot the 10th eigenvalue
plt.plot(domain_f, [np.sin(np.pi*10*x) for x in domain_f])
plt.plot(domain, [np.sin(np.pi*10*x) for x in domain], "+", markersize=20)
plt.show()