import matplotlib.pyplot as plt
import time
from foll import *
from fastar import *
from RHS import *
#from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import numpy as np
from simple_pid import PID

#import warnings
#warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

start_time = time.time()

# ODE solver parameters
abserr = 1.0e-6
relerr = 1.0e-3
stoptime = 100
numpoints = stoptime*10
sec_per_mini_bil = stoptime/numpoints

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.


#t_points = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
t_points = np.linspace(0, stoptime, numpoints+1)


TOW_hradi = 1.25

kp = 1.0
ki = 0.1
kd = 0.0
N = 1.0  # filter coefficient

#controller = PID(setpoint=5*np.pi/180, t_now=37, Kp=2, Ki=0, Kd=0.0)
max_vaeng_horn = 16 * np.pi/180
ref_roll = 5 * np.pi/180
controller = PID(kp, ki, kd, setpoint=ref_roll, sample_time=None,
                 output_limits=(-max_vaeng_horn, max_vaeng_horn))
boo = foo(t_now=37)
x0 = np.array([0,  # u0
               0,  # v0
               0,  # w0
               0,  # p0
               0,  # q0
               0,  # r0
               0,  # x0
               0,  # y0
               0,  # z0
               0,  # phi0
               0,  # theta
               0  # psi0
               ]).reshape(12, 1)  # error integration

w0 = x0.flatten()

# Call the ODE solver.
#wsol = solve_ivp(vectorfield,(0,stoptime), w0, t_eval = t_points, method='BDF', args = (kp, ki, kd, TOW_hradi,controller), atol=abserr, rtol=relerr)
SOL = np.zeros((12,numpoints+1))


L_intv = 2
sec_per_stora_bil = sec_per_mini_bil*(L_intv-1)
styrimerki = 0

for i in range(int((len(t_points)-1)/L_intv)):
    a = L_intv*i
    b = L_intv*i+L_intv
    print(w0)
    
    wsol = solve_ivp(vectorfield, (t_points[a:b][0], t_points[a:b][-1]), w0,
                     t_eval=t_points[a:b], method='BDF', args=(kp, ki, kd, TOW_hradi, styrimerki))

    styrimerki = controller(wsol.y[9,-1],dt=sec_per_stora_bil)
    #print(wsol.y[:])
    #print(wsol.y[:,-1])


    w0 = wsol.y[:,-1]
    SOL[:,a:b] = wsol.y[:]

print('\n')
print("--- execution time: %s seconds ---" % (time.time() - start_time))
print('\n')


# print(wsol)

# print('\n')

# print(wsol.y[0])

plt.figure(1, figsize=(6, 4.5))

plt.xlabel('t')
plt.grid(True)
lw = 1


#plt.plot(wsol.t, wsol.y[9]*180/np.pi, 'b', linewidth=lw)

plt.plot(t_points, SOL[9,:]*180/np.pi, 'b', linewidth=lw)


# print(wsol.y[9]*180/np.pi)


""" fig, axs = plt.subplots(4, 3)
axs[0, 0].plot(t, wsol[:, 0])
axs[0, 0].set_title("u")
axs[0, 1].plot(t, wsol[:, 1])
axs[0, 1].set_title("v")
axs[0, 2].plot(t, wsol[:, 2])
axs[0, 2].set_title("w")

axs[1, 0].plot(t, wsol[:, 3]*180/np.pi)
axs[1, 0].set_title("L_intv")
axs[1, 0].sharex(axs[0, 0])
axs[1, 1].plot(t, wsol[:, 4]*180/np.pi)
axs[1, 1].set_title("q")
axs[1, 1].sharex(axs[0, 1])
axs[1, 2].plot(t, wsol[:, 5]*180/np.pi)
axs[1, 2].set_title("r")
axs[1, 2].sharex(axs[0, 2])

axs[2, 0].plot(t, wsol[:, 6])
axs[2, 0].set_title("x")
axs[2, 0].sharex(axs[0, 0])
axs[2, 1].plot(t, wsol[:, 7])
axs[2, 1].set_title("y")
axs[2, 1].sharex(axs[0, 1])
axs[2, 2].plot(t, wsol[:, 8])
axs[2, 2].set_title("z")
axs[2, 2].sharex(axs[0, 2])

axs[3, 0].plot(t, wsol[:, 9]*180/np.pi)
axs[3, 0].set_title("phi")
axs[3, 0].sharex(axs[0, 0])
axs[3, 1].plot(t, wsol[:, 10]*180/np.pi)
axs[3, 1].set_title("theta")
axs[3, 1].sharex(axs[0, 1])
axs[3, 2].plot(t, wsol[:, 11]*180/np.pi)
axs[3, 2].set_title("psi")
axs[3, 2].sharex(axs[0, 2])

fig.tight_layout() """

plt.show()
