import matplotlib.pyplot as plt
from foll import *
from fastar import *
from RHS import *
import numpy as np
import control as ct

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

S = np.load('stong.npy')
R = np.load('rammi.npy')
t = np.load('timi.npy')

plt.figure()
plt.legend(['phi'])
plt.plot(t, S[9]*180/np.pi,'r')
plt.plot(t, R[9]*180/np.pi,'b')
plt.legend()
plt.grid(True)
plt.show()