from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from ctypes import *
import time
import matplotlib.pyplot as plt
print('subSystem 1')
#
# EMPIRE_API & connect
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
libempire_api.EMPIRE_API_Connect("subSystem1.xml")
#
# This is a linear mass-spring system integrated with BE
# A Dirichlet Neumann coupling of spring vs spring example
# Please see: "Interface-Jacobian based Co-Simulation" Sicklinger et al.
#             International Journal for Numerical Methods in Engineering
# Problem parameters
# Time Step
h = 0.01
# Mass
m = 0.45
# Stiffness
k = 0.5
# Time information
timeSteps = 50
deltaT = 0.1
#
interfaceForceU = (c_double * 1)(0)
interfaceDisplacementY = (c_double * 1)(0)
size = c_int(1)
signal2 = (c_double * 1)(10000)
size2 = c_int(1)
isConvergent = c_int(1)
isConvergent = 0
# Initial condition v0=1 u0=0
u_n_m_1 = [-1 * deltaT]
u_n = [0]
u_n_p_1 = [0]
u_history = [0]
# Iteration counter
iter = 1
# Input-Ouput history
io_history = [0]
#


def doSolve(input, m, k, deltaT, u_n_p_1, u_n, u_n_m_1):
    fac = 1 / (m + deltaT * deltaT * k)
    u_n_p_1[0] = deltaT * deltaT * fac * input + m*fac*(2*u_n-u_n_m_1)
    return u_n_p_1[0]
#


def goToNextTimeStep(u_n_p_1, u_n, u_n_m_1, u_history):
    u_history.append(u_n_p_1[0])
    u_n_m_1[0] = u_n[0]
    u_n[0] = u_n_p_1[0]
    return
#
for i in range(1, timeSteps + 1):
    print('################## Current time step: {}'.format(i))
    while True:
        print('--------------Iteration: {}'.format(iter))
        print('Receiving ...')
        libempire_api.EMPIRE_API_recvSignal_double("interfaceForceU", size, interfaceForceU)
        print('Received interfaceForceU: {}'.format(interfaceForceU[0]))
        io_history.append(interfaceForceU[0])

        interfaceDisplacementY[0] = doSolve(interfaceForceU[0], m, k, deltaT, u_n_p_1, u_n[0], u_n_m_1[0])

        print('Sending ...')
        libempire_api.EMPIRE_API_sendSignal_double("interfaceDisplacementY", size, interfaceDisplacementY)
        print('Sent interfaceDisplacementY: {}'.format(interfaceDisplacementY[0]))
        io_history.append(interfaceDisplacementY[0])

        isConvergent = libempire_api.EMPIRE_API_recvConvergenceSignal()
        print('isConvergent: {}'.format(isConvergent))
        if isConvergent != 0:
            break
        iter = iter + 1
    goToNextTimeStep(u_n_p_1, u_n, u_n_m_1, u_history);
    iter = 1;
    print('u_n_p_1: %f' % u_n_p_1[0])
    print('u_n: %f' % u_n[0])
    print('u_n_m_1: %f' % u_n_m_1[0])

#
# EMPIRE disconnect
libempire_api.EMPIRE_API_Disconnect()
#
# Plot the result
plt.plot(u_history)
# plt.show()
exit();
