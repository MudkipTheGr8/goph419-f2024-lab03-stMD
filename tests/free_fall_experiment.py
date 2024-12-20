# Question 1: Stimulating free fall experiment for heights 10m, 20m, and 40m
#
#
#
#Paeameters:
g0 = 9.811636
gd_dz = 3.086e-6
heights = [10, 20, 40]
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
plt.show()