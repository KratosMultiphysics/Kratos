from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from ctypes import *
import time
print('subSystem 2')
#
# EMPIRE_API & connect
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
libempire_api.EMPIRE_API_Connect("subSystem2.xml")
#
# This is a linear mass-spring system integrated with BE
# A Dirichlet Neumann coupling of spring vs spring example
# Please see: "Interface-Jacobian based Co-Simulation" Sicklinger et al.
#             International Journal for Numerical Methods in Engineering
# Problem parameters
# Time Step
h = 0.01
# Mass
m = 0.55
# Stiffness
k = 0.5
timeSteps = 50
deltaT = 0.1
#

interfaceDisplacementU = (c_double * 1)(0)
interfaceForceY = (c_double * 1)(0)
size = c_int(1)
isConvergent = c_int(1)
isConvergent = 0
# Initial condition v0=1 u0=0
u_n_m_1 = [-1 * deltaT]
u_n = [0]
u_n_p_1 = [0]
# Iteration counter
iter = 1
# Input-Ouput history
io_history = [0]
#


def doSolve(input, m, k, deltaT, u_n_p_1, u_n, u_n_m_1):
    fac = 1 / (m + deltaT * deltaT * k)
    u_n_p_1[0] = input
    f_n_p_1 = -(u_n_p_1[0] * ((m / (deltaT * deltaT)) + k) - (m / (deltaT * deltaT)) * (2 * u_n - u_n_m_1))
    return f_n_p_1
#


def goToNextTimeStep(u_n_p_1, u_n, u_n_m_1):
    u_n_m_1[0] = u_n[0]
    u_n[0] = u_n_p_1[0]
    return
#

for i in range(1, timeSteps + 1):
    print('################## Current time step: {}'.format(i))
    while True:
        print('--------------Iteration: {}'.format(iter))
        print('Receiving ...')
        libempire_api.EMPIRE_API_recvSignal_double("interfaceDisplacementU", size, interfaceDisplacementU)
        print('Received DisplacementU: {}'.format(interfaceDisplacementU[0]))
        io_history.append(interfaceDisplacementU[0])

        interfaceForceY[0] = doSolve(interfaceDisplacementU[0], m, k, deltaT, u_n_p_1, u_n[0], u_n_m_1[0])

        print('Sending ...')
        libempire_api.EMPIRE_API_sendSignal_double("interfaceForceY", size, interfaceForceY)
        print('Sent interfaceForceY: {}'.format(interfaceForceY[0]))
        io_history.append(interfaceForceY[0])

        isConvergent = libempire_api.EMPIRE_API_recvConvergenceSignal()
        print('isConvergent: {}'.format(isConvergent))
        if isConvergent != 0:
            break
        iter = iter + 1
    goToNextTimeStep(u_n_p_1, u_n, u_n_m_1);
    iter = 1;
    print('u_n_p_1: %f' % u_n_p_1[0])
    print('u_n: %f' % u_n[0])
    print('u_n_m_1: %f' % u_n_m_1[0])
#
# Generate data for compare
f = open('io_history_subSystem2.dat', 'w+')
for s in io_history:
    f.write('%e \n' % s)
f.close()
#
# EMPIRE disconnect
libempire_api.EMPIRE_API_Disconnect()
exit();
