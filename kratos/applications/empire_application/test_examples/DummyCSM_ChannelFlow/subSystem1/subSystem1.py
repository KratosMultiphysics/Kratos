import os
from ctypes import *
print 'subSystem 1'

############################################################################
# EMPIRE_API & connect
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
libempire_api.EMPIRE_API_Connect("subSystem1.xml")


#---------------------------------------------------------------------------
interfaceDisplacementY = (c_double * 7)(0)
interfaceForceU	       = (c_double * 7)(0)
size                   = 7
isConvergent  	       = 0

# Time information
timeSteps=10
deltaT=0.1

#loop information
it = 0
maxit = 50

# prescribed rigid-body-motion:
interfaceDisplacementY = [ 0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 ]
interfaceDisplacementY = (c_double * len(interfaceDisplacementY))(*interfaceDisplacementY)

for i in range(1,timeSteps+1):
	print '################## Current time step: {}'.format(i)
	while it < maxit:
       		print '--------------Iteration: {}'.format(it)
		print 'Receiving ...'

		libempire_api.EMPIRE_API_recvSignal_double("interfaceForceU", size, interfaceForceU)
		print 'Received interfaceForceU: {}'.format(interfaceForceU[0])

		print 'Sending ...'
		
		libempire_api.EMPIRE_API_sendSignal_double("interfaceDisplacementY", size, interfaceDisplacementY)
		print 'Sent interfaceDisplacementY: {}'.format(interfaceDisplacementY[3])

       	  	isConvergent=libempire_api.EMPIRE_API_recvConvergenceSignal()
       		print 'isConvergent: {}'.format(isConvergent)
         	if isConvergent!=0:
          		break
		it = it + 1
	it = 1

#---------------------------------------------------------------------------


############################################################################
# EMPIRE disconnect
libempire_api.EMPIRE_API_Disconnect()
exit();

