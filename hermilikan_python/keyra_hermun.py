import matplotlib.pyplot as plt
from matplotlib import animation
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



heildartimi = 6*60 # sek
number_of_solution_points = heildartimi*30

X0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])               # Initial H, L
T = np.linspace(0, heildartimi, number_of_solution_points)   # Simulation time

# PID controller parameters
kp = 1.0
ki = 0.0
kd = 0.0
N = 1.0  # filter coefficient


# function that returns list of multiple step functions
def step_functions(t, t_step, y_step):
    """roll_ref sem fall af tima

    Args:
        t (list): time list
        t_step (list): timi thar sem step kickar inn
        y_step (list): staerd steps

    Returns:
        list: compound step function
    """
    step_functions = []
    for i in range(len(t_step)):
        step_functions.append(np.piecewise(t, [t < t_step[i], t >= t_step[i]], [0, y_step[i]]))
    # sum all lists in list
    step_sum = np.sum(step_functions, axis=0)
    return step_sum

# function that plots a list with respect to time
def plot_time(t, list, title):
    plt.plot(t, list)
    plt.title(title)
    plt.grid(True)
    plt.show()

# Simulate the system
temp = True
if temp == True:
    t1 = 37
    h1 =5./180. * np.pi

    t2 = 149
    h2 =5./180. * np.pi

    t3 = 190
    h3 =5./180. * np.pi

    t4 = 236
    h4 =-15./180. * np.pi

    t_step = [t1,t2,t3,t4]
    y_step = [h1,h2,h3,h4]
    # refrence values fyrir roll (phi) a vorpunni
    roll_ref = step_functions(T, t_step, y_step)
else:
    t1 = 60
    h1 =5./180. * np.pi
    #roll_ref=[0 if t <= t1 else h1 * (t-t1) if t <= t1+1  else h1 for t in T]
    roll_ref = step_functions(T, [t1], [h1])

cl = False

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
    plt.show()
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

#np.save('rammi', y)
#np.save('timi', t)
#b = np.loadtxt('test1.txt', dtype=int)

# function that converts specific indexes of list to degrees
def deg(list, indexes):
    for i in indexes:
        list[i] = list[i] * 180 / np.pi
    return list

y = deg(y, [3,4,5,9,10,11])

# function that plots all the states as a function of time with subplots 3x4
def plot_states(t, y, states):
    fig, axes = plt.subplots(4, 3, figsize=(12, 8))
    for i in range(len(states)):
        axes[i//3, i%3].plot(t, y[i])
        axes[i//3, i%3].set_title(states[i])
        axes[i//3, i%3].grid(True)
        axes[i//3, i%3].sharex(axes[0, i%3])
    fig.tight_layout()
    plt.show()

plot_states(t, y, states_A[0:12])
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')


ax.set_title('3D animation')
#ax.view_init(elev=30, azim=30)
             
def animate(i):
    ax.clear()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim3d(-10+TOW_hradi*t[i], 10+TOW_hradi*t[i])
    ax.set_ylim3d(-10, 10)
    ax.set_zlim3d(-10, 10)
    #ax.clear()
    [J, R, T_tf] = eulerang(y[9][i], y[10][i], y[11][i])
    o_b = np.array([y[6][i],y[7][i],-y[8][i]]).reshape(3,1)
    A = (o_b+R@(R_HO)).flatten()
    B = (o_b+R@(R_HN)).flatten()
    C = (o_b+R@(R_VO)).flatten()
    D = (o_b+R@(R_VN)).flatten()
    #print(A)
    #print(B)
    #print(C)
    #print(D)
    ax.plot([A[0],B[0]], [A[1],B[1]], [A[2],B[2]], 'b')
    ax.plot([A[0],C[0]], [A[1],C[1]], [A[2],C[2]], 'r')
    ax.plot([C[0],D[0]], [C[1],D[1]], [C[2],D[2]], 'g')
    ax.plot([B[0],D[0]], [B[1],D[1]], [B[2],D[2]], 'r')
    

    
    #plt.pause(0.01)
    #plt.autoscale(enable=True, axis='both')

# function that uses matplotlib.animation to animate the lists with respect to time

R_HO[2,0] = -R_HO[2,0]
R_HN[2,0] = -R_HN[2,0]
R_VO[2,0] = -R_VO[2,0]
R_VN[2,0] = -R_VN[2,0]

ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=1)
plt.show()