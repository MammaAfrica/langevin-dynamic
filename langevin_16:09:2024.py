import numpy as np
import matplotlib.pyplot as plt
import math as mt

# the equation that is used in this code is : v(t)=sqrt(2)*eta(t) where eta(t) is a wiener process with 
# square of the standard deviation equal to the timestep of integration
# the diffusion coefficient that is obtained should be one for every dimension

# selection of parameters
D=1                   # number of system dimensions
N=1000                # number of particles in the simulation box
time_step=0.01        # timestep implemented in the integration
num_iter=10000        # numero di iterazioni

#functions
def integration_position(N,D,r,v,time_step):
    for i in range(N):
        for j in range(D):
            r[i,j]=r[i,j]+v[i,j]*time_step
    return r
    
def integration_velocity(N,D,v,friction,time_step,Gamma, mass):
    for i in range(N):
        for j in range(D):
            v[i,j]=v[i,j]-friction*v[i,j]*time_step/mass
            v[i,j]=v[i,j]+Gamma*np.random.normal(0, mt.sqrt(2*time_step), 1)[0]/mass
    return v

def calc_MSD(N,D,r):
    r2=0
    for i in range(N):  # updating the mean square displacement 
        for j in range(D):
           r2=r2+r[i,j]**2
    r2=r2/N
    return r2

def inizialization(D,N,num_iter):
    # inizialization of velocities and positions
    r=np.zeros((N, D))
    v=np.zeros((N, D))
    r2=np.zeros(num_iter)
    return r,v,r2

def integration_loop(D, N, num_iter, time_step, friction, Gamma, mass):
    r,v,r2 =inizialization(D,N,num_iter)
    for i in range(num_iter):
        r=integration_position(N,D,r,v,time_step)
        v=integration_velocity(N,D,v,friction,time_step,Gamma, mass)
        r2[i]= calc_MSD(N,D,r)
    return r2


mass=np.array([2, 1.5, 1, 0.5])
friction=1
Gamma=1

# plot
fig, ax=plt.subplots()
time=np.arange(num_iter)*time_step

for i in range(3):

    r2=integration_loop(D, N, num_iter, time_step, friction, Gamma, mass[i])
    ax.plot(time, r2)


ax.set_xlabel('Time')
ax.set_ylabel('MSD/2/D')
ax.legend(mass)
plt.show()