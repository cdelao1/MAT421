import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
import scipy.special as sp1
from mpl_toolkits.mplot3d import Axes3D
import time


#Project Part 1
start = time.time()

def v(r,z,gamma):
    a=r*(1-z/gamma) 
    sums = 0
    for n in range(1,int(gamma*200)):
        sums += ((sp1.iv(1,(n*np.pi*r)/gamma))/(n*sp1.iv(1,(n*np.pi)/gamma)))*np.sin(n*np.pi*z/gamma)
    return a-(2/np.pi)*sums
    
def plot_contour(a, filename=None, zlabel='v(r,z)',cmap=plt.cm.gnuplot):
    fig = plt.figure(figsize=(4.5,4.5))
    ax = fig.add_subplot(111)

    x = np.arange(a.shape[0])
    y = np.arange(a.shape[1])
    X, Y = np.meshgrid(x, y)
    Z = a[X, Y]
    cset = ax.contourf(X, Y, Z, 20, cmap=cmap)
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_title('\u0393=1.5')
    ax.axis('off')
    ax.set_aspect(1)

    cb = fig.colorbar(cset, shrink=0.5, aspect=5)
    cb.set_label(zlabel)
    
    if filename:
        fig.savefig(filename,dpi=300)
        plt.close(fig)
        return filename
    else:
        return ax

gamma = 1.5
Nmax = 100
M = int(Nmax*gamma)

v1 = np.zeros((Nmax+2,M+2),dtype=np.float64)
for r in range(Nmax+2):
    for z in range(M+2):
        v1[r][z] = v(1/(Nmax+1)*r,1/(M+1)*z*gamma,gamma)

plot_contour(v1, 'gamma15e')
print(time.time()-start)

##########################################################
#Project Part 2

import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
import scipy.special as sp1
from mpl_toolkits.mplot3d import Axes3D
import time

start = time.time()

def v(r,z,gamma):
    a=r*(1-z/gamma) 
    sums = 0
    for n in range(1,int(gamma*200)):
        sums += ((sp1.iv(1,(n*np.pi*r)/gamma))/(n*sp1.iv(1,(n*np.pi)/gamma)))*np.sin(n*np.pi*z/gamma)
    return a-(2/np.pi)*sums
    
def plot_contour(a, filename=None, zlabel='v(r,z)',cmap=plt.cm.gnuplot):
    fig = plt.figure(figsize=(4.5,4.5))
    ax = fig.add_subplot(111)

    x = np.arange(a.shape[0])
    y = np.arange(a.shape[1])
    X, Y = np.meshgrid(x, y)
    Z = a[X, Y]
    cset = ax.contourf(X, Y, Z, 20, cmap=cmap)
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_title('\u0393=1.5')
    ax.axis('off')
    ax.set_aspect(1)

    cb = fig.colorbar(cset, shrink=0.5, aspect=5)
    cb.set_label(zlabel)
    
    if filename:
        fig.savefig(filename,dpi=300)
        plt.close(fig)
        return filename
    else:
        return ax

gamma = 1.5
Nmax = 100
M = int(Nmax*gamma)
        
A_nn = np.zeros((Nmax, Nmax), dtype=np.float64)
#main diagonal
for i in range(Nmax):
    dr = 1/(Nmax)
    A_nn[i,i] = -(2+(1/((i+1)**2)))/(dr**2)
#sub diagonal
for i in range(Nmax-1):
    dr = 1/(Nmax)
    A_nn[i+1,i] = (1-1/(2*(i+1)))/(dr**2)
#super diagonal
for i in range(Nmax-1):
    dr = 1/(Nmax)
    A_nn[i,i+1] = (1+1/(2*(i+1)))/(dr**2)
    
e_i, Z_nn = np.linalg.eig(A_nn)
invZ_nn = np.linalg.inv(Z_nn)

F_nm = np.zeros((Nmax, M), dtype=np.float64)
for n in range(Nmax):
    dz = gamma/M
    dr = 1/Nmax  
    F_nm[n,0] = -(n+1)*(dr/(dz**2))
 
B_mm = np.zeros((M,M), dtype = np.float64)
#main diagonal
for j in range(M):
    dz = gamma/M
    B_mm[j,j] = -2/(dz**2)
#sub diagonal
for j in range(M-1):
    dz = gamma/M
    B_mm[j+1,j] = 1/(dz**2)
#super diagonal
for j in range(M-1):
    dz = gamma/M
    B_mm[j,j+1] = 1/(dz**2)

H_mn = np.transpose(F_nm).dot(np.transpose(invZ_nn))
I_mm = np.identity(M)
U_nm = np.zeros((Nmax, M), dtype=np.float64)
for i in range(Nmax):
    U_nm[i,:] = np.linalg.inv(B_mm + e_i[i]*I_mm).dot(H_mn[:,i])

v2 = Z_nn.dot(U_nm)
plot_contour(v2, 'gamma15n')

print(time.time()-start)

##########################################################
#Project Part 3
import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
import scipy.special as sp1
from mpl_toolkits.mplot3d import Axes3D
import time

start = time.time()

def plot_contour(a, filename=None, zlabel='v(r,z)',cmap=plt.cm.gnuplot):
    fig = plt.figure(figsize=(4.5,4.5))
    ax = fig.add_subplot(111)
    ax 

    x = np.arange(a.shape[0])
    y = np.arange(a.shape[1])
    X, Y = np.meshgrid(x, y)
    Z = a[X, Y]
    cset = ax.contourf(X, Y, Z, 20, cmap=cmap)
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_title('\u0393=1.5')
    ax.axis('off')
    ax.set_aspect(1)

    cb = fig.colorbar(cset, shrink=0.5, aspect=5)
    cb.set_label(zlabel)
    
    fig.savefig(filename,dpi=300)
    plt.close(fig)
    return ax

gamma = 1.5
n = 50
m = int(gamma*n)
t=0
dz = 1/n
dr = 1/n
dt = 1e-4
Re = 1
tol = 1e-5
stop = 15000
r = np.linspace(0,1,n)

def L2(v, nr=n, nz=n):
    v_norm = np.sqrt(1/((nr+1)*(nz+1))*np.sum(np.power(v,2)))
    return v_norm

def RHS(v):
    u = np.zeros((n,m))
    for i in range(1,n-1):
        for j in range(1,m-1):
            a = (v[i+1,j]-2*v[i,j]+v[i-1,j])
            b = (v[i+1,j]-v[i-1,j])/(2*i)
            c = v[i,j]/i**2
            d = (v[i,j+1]-2*v[i,j]+v[i,j-1])
            u[i,j] = (1/Re)*(a+b-c+d)/dr**2
    return(u)

def Heun(v):
    pre = v+dt*RHS(v)
    return v+0.5*dt*(RHS(v)+RHS(pre))

def limit(v1, v2, lim = tol):
    if (L2(v1)-L2(v2))/L2(v1) < lim:
        return True
    else:
        return False

v0 = np.zeros((n,m))
v0[:,0] = r
v1 = Heun(v0)


while limit(v0,v1) == True:
    v0 = Heun(v0)
    v1 = Heun(v1)
    if t < 100:
        if t % 5 == 0:
            plot_contour(v0, filename=str(t))
    elif t < 250:
        if t % 10 == 0:
            plot_contour(v0, filename=str(t))
    elif t < 1000:
        if t % 50 == 0:
            plot_contour(v0, filename=str(t))
    else:
        if t % 100 == 0:
            plot_contour(v0, filename=str(t))
    t += 1

y = time.time()-start
print("This took", "{:.2f}".format(y/60), "minutes to run!")