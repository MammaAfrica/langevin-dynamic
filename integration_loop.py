import numpy as np
import matplotlib.pyplot as plt
from funzioni.corr_MSD import MSD
from funzioni.corr_A_B import corr_A_B
from funzioni.inizialization import inizialization
from funzioni.update_position import update_position
from funzioni.update_velocity import update_velocity
from funzioni.update_velocity import force_calc



time_step=0.01         # timestep implemented in the integration
num_iter=50000       # numero di iterazioni
num_corr=100
num_corr_MSD=500

sigma=1
L_box=20*sigma


mass=0.1
friction=1
Gamma=1

r0=5
v0=0
potential_type="sin"
potential="yes"
boundary="periodic"


r, v, corr_data, corr_avg, n_corr_points, corr_data_MSD, corr_avg_MSD, n_corr_points_MSD, flag =inizialization(
    r0, v0, num_corr, num_corr_MSD)

for i in range(num_iter):
    r, flag=update_position(r,v,time_step, flag, boundary, L_box)

    v, noise, Force=update_velocity(r,v, friction ,time_step, Gamma, mass, 
                                    potential_type, potential, sigma, L_box )
    
    # <B(0)A(-t)>
    B=Force
    A=noise
    
    corr_data, corr_avg, n_corr_points = corr_A_B(corr_data.copy(), corr_avg.copy(),
                                                   n_corr_points.copy(), num_corr, A, B, i)
    
    corr_data_MSD, corr_avg_MSD, n_corr_points_MSD=MSD(corr_data_MSD.copy(), 
                                                       corr_avg_MSD.copy(), n_corr_points_MSD.copy(), 
                                                       num_corr_MSD, r, flag, i, L_box)
    
    if i%10000==0:
        print(flag)

#plot

delay_MSD=np.arange(num_corr_MSD)*time_step
delay=np.arange(num_corr)*time_step


fig1, ax=plt.subplots(1)
ax.plot(delay_MSD, corr_avg_MSD/2)
ax.plot(delay_MSD, delay_MSD)
ax.set_xlabel('time')
ax.set_ylabel('MSD')

fig1, ax1=plt.subplots(1)

ax1.plot(delay, corr_avg)
ax1.set_xlabel('time')
ax1.set_ylabel('<v(0)v(t)>')

plt.show()