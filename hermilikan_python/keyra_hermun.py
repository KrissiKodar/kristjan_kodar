import matplotlib.pyplot as plt
import time
from foll import *
from fastar import *
from RHS import *
from scipy.integrate import odeint
import numpy as np
#import warnings
#warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

start_time = time.time()

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 300
numpoints = 10000

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.


t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
w0 = x0.flatten().tolist()

# Call the ODE solver.
wsol = odeint(vectorfield, w0, t, args = (kp, ki, kd,), atol=abserr, rtol=relerr)

print('\n')
print("--- execution time: %s seconds ---" % (time.time() - start_time))
print('\n')
plt.figure(1, figsize=(6, 4.5))

plt.xlabel('t')
plt.grid(True)
lw = 1

plt.plot(t, wsol[:, 9]*180/np.pi, 'b', linewidth=lw)



fig, axs = plt.subplots(4, 3)
axs[0, 0].plot(t, wsol[:, 0])
axs[0, 0].set_title("u")
axs[0, 1].plot(t, wsol[:, 1])
axs[0, 1].set_title("v")
axs[0, 2].plot(t, wsol[:, 2])
axs[0, 2].set_title("w")

axs[1, 0].plot(t, wsol[:, 3]*180/np.pi)
axs[1, 0].set_title("p")
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

fig.tight_layout()

plt.show()



