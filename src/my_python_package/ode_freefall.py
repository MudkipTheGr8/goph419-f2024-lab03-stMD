# Name: Matthew Davidson UCID: 30182729
#
#Purpose:Solving second-order ODE through using Euler's Method and classical 4th order RungeKutta (RK4)
#
#Variables:
# g0: gravitational acceleration
# dg_dz: free air gradient of gravitational acceleration
# cd_star: mass normalized drag coefficient
# H: Total drop height
# dt: time step for the ODE solver
# t: time
# z: distance travelled
# v: velocity
#
#
import time as tm
import numpy as np
import matplotlib.pyplot as plt
#
#Parameters: g0, dg_dz, cd_star, H, dt
#Returns: t, z, v
#
def ode_freefall_euler(g0, dg_dz, cd_star, H, dt):
    # Initialization of time, height, and velocity
    t = [0]
    z = [0]
    v = [0]
    # A loop will occur until the height is reached 
    while z[-1] < H:
        new_z = z[-1]
        new_v = v[-1]
        #calculating acceleration
        acc = g0 - cd_star * new_v - dg_dz * new_z
        #Euler's Method
        v_up = new_v + acc * dt
        z_up = new_z + new_v * dt
        #append new position, velocity, and time
        z.append(z_up)
        v.append(v_up)
        t.append(t[-1] + dt)
        # Confirming position equals H
        if z[-1] > H:
            extra = z[-1] - H
            dt_ad = dt - (extra / new_v)
            z[-1] = H
            t[-1] = t[-2] + dt_ad
    
    return np.array(t), np.array(z), np.array(v)

#Parameters: g0, dg_dz, cd_star, H, dt
#Returns: t, z, v
#
def ode_freefall_rk4(g0, dg_dz, cd_star, H, dt):
    # Initialization of time, height, and velocit
    t = [0]
    z = [0]
    v = [0]
    # A loop will occur until the height is reached
    while z[-1] < H:
        new_z = z[-1]
        new_v = v[-1]
        #calculating k and l for RungeKutta
        k1 = new_v
        l1 = g0 - cd_star * new_v - dg_dz * new_z
        k2 = new_v + (1 / 2) * dt * l1
        l2 = g0 - cd_star * (new_v + (1/2) * dt * l1) - dg_dz * (new_z + (1/2) * dt * k1)
        k3 = new_v + (1/2) * dt * l2
        l3 = g0 - cd_star * (new_v + (1/2) * dt * l2) - dg_dz * (new_z + (1/2) * dt * k2)
        k4 = new_v + dt * l3
        l4 = g0 - cd_star * (new_v + (1/2) * dt * l3) - dg_dz * (new_z + (1/2) * dt * k3)
        #RK4 Method
        z_up = new_z + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        v_up = new_v + (dt / 6) * (l1 + 2 * l2 + 2 * l3 +l4)
        # Appending into the list
        z.append(z_up)
        v.append(v_up)
        t.append(t[-1] + dt)
        #Confirming position equals H
        if z[-1] > H:
            extra = z[-1] - H
            dt_ad = dt - (extra / new_v)
            z[-1] = H
            t[-1] = t[-2] + dt_ad
    return np.array(t), np.array(z), np.array(v)
        



# Question 1: Stimulating free fall experiment for heights 10m, 20m, and 40m
#
#
#
#Paeameters:
g0 = 9.811636
dg_dz = 3.086e-6
cd_star = 0.5
heights = [10, 20, 40]
dt = 0.01
#
#time steps:
time = np.linspace(0.001, 0.1, 20)
#
#Storing results:
#
euler_time = {H: [] for H in heights}
RK4_time = {H: [] for H in heights}
euler_error = {H: [] for H in heights}
RK4_error = {H: [] for H in heights}
#
#Stimulating free fall
#
for H in heights:
    for dt in time:
        #Euler Method
        euler_t, euler_z, euler_v = ode_freefall_euler(g0, dg_dz, cd_star, H, dt)
        euler_time[H].append(euler_t[-1])
        euler_error[H].append(abs((euler_t[-1] - euler_t[-2])/euler_t[-1]))
        #RK4 Method
        rk4_t, rk4_z, rk4_v = ode_freefall_rk4(g0, dg_dz, cd_star, H, dt)
        RK4_time[H].append(rk4_t[-1])
        RK4_error[H].append(abs((rk4_t[-1] - rk4_t[-2]) / rk4_t[-1]))

# forming graph
plt.figure(figsize=(14,6))
plt.subplot(1, 2, 1)
for H in heights:
    plt.plot(time, euler_time[H], label=f'Euler H={H}m')
    plt.plot(time, RK4_time[H],'--', label=f'RK4 H={H}m')
plt.xlabel('Time step(s)')
plt.ylabel('Total drop time(s)')
plt.title('Total Drop Time vs Time Step')
plt.legend()
plt.grid(True)

# forming error graph
plt.subplot(1, 2, 2)
for H in heights:
    plt.plot(time, euler_error[H], label=f'Euler H={H}m')
    plt.plot(time, RK4_error[H], '--', label=f'RK4 H={H}m')
plt.xlabel('Time step(s)')
plt.ylabel('Approximate Relative Error')
plt.title('Relative Error vs Time Step')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('/Users/matth/OneDrive/Desktop/GOPH_LAB_3/goph419-f2024-lab03-stMD/figures/Question_1_Graphs.png',bbox_inches="tight")


# Question #2
#
#

euler_time = {H: [] for H in heights}
rk4_time = {H: [] for H in heights}
for H in heights:
    for dt in time:
        #Euler Method
        start = tm.perf_counter()
        ode_freefall_euler(g0, dg_dz, cd_star, H, dt)
        end = tm.perf_counter()
        euler_time[H].append(end - start)
        # RK4 Method
        start = tm.perf_counter()
        ode_freefall_rk4(g0, dg_dz, cd_star, H, dt)
        end = tm.perf_counter()
        rk4_time[H].append(end - start)

#Euler 
plt.figure(figsize=(14, 6))
plt.subplot(1, 2, 1)
for H in heights:
    plt.plot(time, euler_time[H], label=f'Euler H={H}m')
plt.xlabel('Time Step')
plt.ylabel('Simulation Time (s)')
plt.title('Euler Simulation Time vs Time Step')
plt.legend()
plt.grid(True)

#RK4
plt.subplot(1, 2, 2)
for H in heights:
    plt.plot(time, rk4_time[H], label=f'RK4 H={H}m')
plt.xlabel('Time Step')
plt.ylabel('Simulation Time (s)')
plt.title('RK4 Simulation Time vs Time Step')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('/Users/matth/OneDrive/Desktop/GOPH_LAB_3/goph419-f2024-lab03-stMD/figures/Question_2_Graphs.png',bbox_inches="tight")

