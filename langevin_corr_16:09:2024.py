import numpy as np
import matplotlib.pyplot as plt
import math as mt

# the equation that is used in this code is : v(t)=sqrt(2)*eta(t) where eta(t) is a wiener process with 
# square of the standard deviation equal to the timestep of integration
# the diffusion coefficient that is obtained should be one for every dimension

# selection of parameters
D=1                   # number of system dimensions
N=1                   # number of particles in the simulation box
time_step=0.1         # timestep implemented in the integration
num_iter=100000        # numero di iterazioni
num_corr=100

#functions
def integration_position(r,v,time_step):
    r=r+v*time_step
    return r
    
def integration_velocity(r,v,friction,time_step,Gamma, mass):
    v=v-friction*v*time_step/mass
    v=v+Gamma*np.random.normal(0, mt.sqrt(2*time_step), 1)[0]/mass
    return v

def inizialization(num_iter, num_corr):
    # inizialization of velocities and positions
    r=0
    v=0
    # inizialization of the correlate function
    corr_data=np.zeros(num_corr)
    corr_avg=np.zeros(num_corr)
    n_corr_points=np.zeros(num_corr)
    return r, v, corr_data, corr_avg, n_corr_points

def integration_loop(num_corr, num_iter, time_step, friction, Gamma, mass):

    r, v, corr_data, corr_avg, n_corr_points =inizialization(num_iter, num_corr)

    for i in range(num_iter):
        r=integration_position(r,v,time_step)
        v=integration_velocity(r,v,friction,time_step,Gamma, mass)
        corr_data, corr_avg, n_corr_points = correlate_fly(corr_data.copy(), corr_avg.copy(), n_corr_points.copy(), num_corr, num_iter, v, i)
    return  corr_avg

def correlate_fly(corr_data, corr_avg, n_corr_points, num_corr, num_iter, v, integration_step):
    # update the data vector
    
    corr_data_0=corr_data.copy()
    corr_data[0]=v
    for i in range(1, num_corr):
          corr_data[i]=corr_data_0[i-1]

    # update the correlation functuion
    for i in range(num_corr):
        if integration_step >=i:
            corr_avg[i]=(corr_avg[i]*n_corr_points[i]+corr_data[0]*corr_data[i])/(n_corr_points[i]+1)
            n_corr_points[i]=n_corr_points[i]+1

    return corr_data, corr_avg, n_corr_points


mass=np.array([0.1, 0.5, 1])
friction=1
Gamma=1


#plot
fig, ax=plt.subplots(1)
delay=np.arange(num_corr)*time_step

for i in range(3):
    corr_fun =integration_loop(num_corr, num_iter, time_step, friction, Gamma, mass[i])
    ax.scatter(delay, corr_fun)

ax.set_xlabel('delay')
ax.set_ylabel('corr_fun')

plt.show()