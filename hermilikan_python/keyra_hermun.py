import matplotlib.pyplot as plt
import time
from foll import *
from fastar import *
from RHS import *
from scipy.integrate import solve_ivp
import numpy as np
import control as ct

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

start_time = time.time()


################################ plant #########################################

states  = ('ub','vb','wb','p','q','r','x','y','z','phi','theta','psi')
outputs = ('ub','vb','wb','p','q','r','x','y','z','phi','theta','psi')
io_varpa= ct.NonlinearIOSystem(
    varpa_dynamics, None, inputs=('u'), outputs=outputs,
    states=states, name='varpa')

############################### controller ######################################


def pi_update(t, x, u, params={}):
    # Get the controller parameters that we need
    ki = params.get('ki', 0.1)
    kaw = params.get('kaw', 2)  # anti-windup gain

    # Assign variables for inputs and states (for readability)
    v = u[0]                    # current velocity
    vref = u[1]                 # reference velocity
    z = x[0]                    # integrated error

    # Compute the nominal controller output (needed for anti-windup)
    u_a = pi_output(t, x, u, params)

    # Compute anti-windup compensation (scale by ki to account for structure)
    u_aw = kaw/ki * (np.clip(u_a, 0, 1) - u_a) if ki != 0 else 0

    # State is the integrated error, minus anti-windup compensation
    return (vref - v) + u_aw

def pi_output(t, x, u, params={}):
    # Get the controller parameters that we need
    kp = params.get('kp', 0.5)
    ki = params.get('ki', 0.1)
    kd = params.get('kd', 0.1)
    # Assign variables for inputs and states (for readability)
    v = u[0]                    # current velocity
    vref = u[1]                 # reference velocity
    v_dot = u[2]                # hradi snunings (nokkurn veginn)
    z = x[0]                    # integrated error

    # PI controller
    return kp * (vref - v) + ki * z - kd * v_dot

control_pi = ct.NonlinearIOSystem(
    pi_update, pi_output, name='control',
    inputs = ['phi', 'roll_ref','p'], outputs = ['u'], states = ['z'])

outputs_A = ['ub','vb','wb','p','q','r','x','y','z','phi','theta','psi','u']
states_A = ['ub','vb','wb','p','q','r','x','y','z','phi','theta','psi','u']
# Create the closed loop system for the state space controller

hermun = ct.InterconnectedSystem(
    (io_varpa, control_pi), name='hermun',connections=(('varpa.u', 'control.u'),('control.phi', 'varpa.phi'),('control.p', 'varpa.p')),
        inplist=('control.roll_ref'),
        outlist=('varpa.ub','varpa.vb','varpa.wb','varpa.p','varpa.q','varpa.r','varpa.x','varpa.y','varpa.z','varpa.phi','varpa.theta','varpa.psi','control.u'),
        outputs=outputs_A)


X0 = np.zeros(12)#[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]                 # Initial H, L
T = np.linspace(0, 6*60, 5000)   # Simulation 70 years of time



kp = 1.0
ki = 0.0
kd = 0.0
N = 1.0  # filter coefficient


# Simulate the system
temp = False
if temp == True:
    t1 = 37
    h1 =5./180. * np.pi

    t2 = 149
    h2 =10./180. * np.pi

    t3 = 190
    h3 =15./180. * np.pi

    t4 = 236

    roll_ref = [0 if t <= t1 else h1 * (t-t1) if t <= t1+1  else h1 
    if t <= t2 else h1 * (t-t2) + h1 if t <= t2+1 else h2 
    if t <= t3 else h1 * (t-t3) + h2 if t <= t3+1 else h3  
    if t <= t4 else 0  for t in T]
else:
    t1 = 60
    h1 =5./180. * np.pi
    roll_ref=[0 if t <= t1 else h1 * (t-t1) if t <= t1+1  else h1 for t in T]

cl = True

# closed loop
if cl == True:
    t, y = ct.input_output_response(hermun, T, roll_ref, X0,solve_ivp_method="LSODA",params = {'kp':kp, 'ki':ki, 'kd':kd})
    # Plot the response
    plt.figure()
    plt.plot(t, y[12]*180/np.pi,'b',label = "AOT vaengur")
    plt.plot(t,  [x * 180/np.pi for x in roll_ref],'g',label = "roll_ref")
    #plt.plot(t, roll_ref*180/np.pi,'g',label = "roll_ref")
    plt.legend(['u'])
    plt.plot(t, y[9]*180/np.pi,'r',label = "roll")
    plt.legend()
    plt.grid(True)
    print(np.shape(y))
else:
    # open loop
    vaeng_horn = 5
    t, y = ct.input_output_response(io_varpa, T, vaeng_horn*np.pi/180, X0,solve_ivp_method="LSODA")
    # Plot the response
    plt.figure()
    plt.legend(['u'])
    plt.plot(t, y[9]*180/np.pi,'r',label = f"roll, u = {vaeng_horn} deg")
    plt.legend()
    plt.grid(True)

plt.show()

print('\n')
print("--- execution time: %s seconds ---" % (time.time() - start_time))
print('\n')






#### sma ad profa ####
""" # xyz_hnit = np.vstack((HO_X, HN_X, VO_X, VN_X))

#rammi_xyz_hnit = np.vstack(
#    (o_b + R@R_HO, o_b + R@R_HN,  o_b + R@R_VO, o_b + R@R_VN))

for i in range(len(t)):
    [J, R, T_tf] = eulerang(y[10][i], y[11][i], y[12][i])
    A=y[7:10][i]+R@(R_HO)
    B=y[7:10][i]+R@(R_HN)
    C=y[7:10][i]+R@(R_VO)
    D=y[7:10][i]+R@(R_VN) """



##### plotta oll state variables #####

fig, axs = plt.subplots(4, 3)
axs[0, 0].plot(t, y[0])
axs[0, 0].set_title("u")
axs[0, 1].plot(t, y[1])
axs[0, 1].set_title("v")
axs[0, 2].plot(t, y[2])
axs[0, 2].set_title("w")

axs[1, 0].plot(t, y[3]*180/np.pi)
axs[1, 0].set_title("p")
axs[1, 0].sharex(axs[0, 0])
axs[1, 1].plot(t, y[4]*180/np.pi)
axs[1, 1].set_title("q")
axs[1, 1].sharex(axs[0, 1])
axs[1, 2].plot(t, y[5]*180/np.pi)
axs[1, 2].set_title("r")
axs[1, 2].sharex(axs[0, 2])

axs[2, 0].plot(t, y[6])
axs[2, 0].set_title("x")
axs[2, 0].sharex(axs[0, 0])
axs[2, 1].plot(t, y[7])
axs[2, 1].set_title("y")
axs[2, 1].sharex(axs[0, 1])
axs[2, 2].plot(t, y[8])
axs[2, 2].set_title("z")
axs[2, 2].sharex(axs[0, 2])

axs[3, 0].plot(t, y[9]*180/np.pi)
axs[3, 0].set_title("phi")
axs[3, 0].sharex(axs[0, 0])
axs[3, 1].plot(t, y[10]*180/np.pi)
axs[3, 1].set_title("theta")
axs[3, 1].sharex(axs[0, 1])
axs[3, 2].plot(t, y[11]*180/np.pi)
axs[3, 2].set_title("psi")
axs[3, 2].sharex(axs[0, 2])

fig.tight_layout()

plt.show()
