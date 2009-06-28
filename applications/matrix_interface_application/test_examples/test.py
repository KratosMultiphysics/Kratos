##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

print "aaa"
#importing Kratos main library
from Kratos import *
print "bbbb"
kernel = Kernel()   #defining kernel
print kernel
print "ccc"

#importing applications
import applications_interface
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosMatrixInterfaceApplication import *

# There is a problem here...

# Of course, this is not correct, as we want a new object created...
m = CompressedMatrix()

# Also, it seems that m is not touched, despite successful reading of the matrix
print ReadMatrixMarket("test-general.mtx", m)
print m
