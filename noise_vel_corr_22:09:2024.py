import numpy as np
import matplotlib.pyplot as plt
import math as mt

# the equation that is used in this code is : v(t)=sqrt(2)*eta(t) where eta(t) is a wiener process with 
# square of the standard deviation equal to the timestep of integration
# the diffusion coefficient that is obtained should be one for every dimension

# selection of parameters
D=1                   # number of system dimensions
N=1                   # number of particles in the simulation box
time_step=0.05         # timestep implemented in the integration
num_iter=100000       # numero di iterazioni
num_corr=100
sigma=1
L_box=10*sigma
num_corr_MSD=200

#functions
def integration_position(r,v,time_step, L_box, sigma):
    r=r+v*time_step
    return r
    
def integration_velocity(r,v,friction,time_step,Gamma, mass):

    v=v-friction*v*time_step/mass+(1/(r-1))**2
    noise=Gamma*np.random.normal(0, mt.sqrt(2*time_step), 1)[0]/mass
    v=v+noise
    return v, noise

def inizialization(num_iter, num_corr, num_corr_MSD):
    # inizialization of velocities and positions
    r=2
    v=0
    # inizialization of the correlate function
    corr_data=np.zeros(num_corr)
    corr_data_MSD=np.zeros(num_corr_MSD)
    corr_avg=np.zeros(num_corr)
    corr_avg_MSD=np.zeros(num_corr_MSD)
    n_corr_points=np.zeros(num_corr)
    n_corr_points_MSD=np.zeros(num_corr_MSD)
    return r, v, corr_data, corr_avg, n_corr_points, corr_data_MSD, corr_avg_MSD, n_corr_points_MSD

def integration_loop(num_corr, num_corr_MSD, num_iter, time_step, friction, Gamma, mass, L_box, sigma):

    r, v, corr_data, corr_avg, n_corr_points, corr_data_MSD, corr_avg_MSD, n_corr_points_MSD =inizialization(num_iter, num_corr, num_corr_MSD)

    for i in range(num_iter):
        r=integration_position(r,v,time_step, L_box, sigma)
        v, noise=integration_velocity(r,v,friction,time_step,Gamma, mass)
        corr_data, corr_avg, n_corr_points = correlate_fly(corr_data.copy(), corr_avg.copy(), n_corr_points.copy(), num_corr, v, noise, i)
        corr_data_MSD, corr_avg_MSD, n_corr_points_MSD= MSD(corr_data_MSD.copy(), corr_avg_MSD.copy(), n_corr_points_MSD.copy(), num_corr_MSD, r, i)
    
    return  corr_avg_MSD, corr_avg

def correlate_fly(corr_data, corr_avg, n_corr_points, num_corr, v, noise, integration_step):
    # update the data vector
    A=noise # earlier in time
    B=v
    corr_data_0=corr_data.copy()
    corr_data[0]=A
    for i in range(1, num_corr):
          corr_data[i]=corr_data_0[i-1]

    # update the correlation functuion
    for i in range(num_corr):
        if integration_step >=i:
            corr_avg[i]=(corr_avg[i]*n_corr_points[i]+B*corr_data[i])/(n_corr_points[i]+1)
            n_corr_points[i]=n_corr_points[i]+1
        
    return corr_data, corr_avg, n_corr_points

def MSD(corr_data_MSD, corr_avg_MSD, n_corr_points_MSD, num_corr_MSD, r, integration_step):
    # update the data vector
    
    corr_data_MSD_0=corr_data_MSD.copy()
    corr_data_MSD[0]=r
    for i in range(1, num_corr_MSD):
          corr_data_MSD[i]=corr_data_MSD_0[i-1]

    # update the correlation functuion
    for i in range(1, num_corr_MSD):
        if integration_step >=i:
            delta=corr_data_MSD[0]-corr_data_MSD[i]
            delta=delta**2
            corr_avg_MSD[i]=(corr_avg_MSD[i]*n_corr_points_MSD[i]+delta)/(n_corr_points_MSD[i]+1)
            n_corr_points_MSD[i]=n_corr_points_MSD[i]+1

    return corr_data_MSD, corr_avg_MSD, n_corr_points_MSD


mass=0.1
friction=1
Gamma=1


#plot
fig, ax=plt.subplots(2)
delay_MSD=np.arange(num_corr_MSD)*time_step
delay=np.arange(num_corr)*time_step

corr_fun_MSD, corr_fun_Vel =integration_loop(num_corr, num_corr_MSD, num_iter, time_step, friction, Gamma, mass, L_box, sigma)
ax[0].plot(delay_MSD, corr_fun_MSD)
ax[1].plot(delay, corr_fun_Vel)

ax[0].set_xlabel('time')
ax[0].set_ylabel('MSD')

ax[1].set_xlabel('time')
ax[1].set_ylabel('<v(0)v(t)>')

plt.show()