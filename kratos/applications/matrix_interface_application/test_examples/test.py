##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

from Kratos import *
kernel = Kernel()   #defining kernel

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosMatrixInterfaceApplication import *

# Create an empty CompressedMatrix
m = CompressedMatrix()

# General case

print "==== Testing general case ========================================="

# Read a test matrix
if (ReadMatrixMarket("test-general.mtx", m)):
	print "Matrix read successfully!"
else:
	print "Error reading matrix!"

# Print the contents of m
print m

if (WriteMatrixMarket("test-general-out.mtx", m, False)):
	print "Matrix written successfully!"
else:
	print "Error writing matrix!"

# Symmetric case

print "==== Testing symmetric case ======================================="

# Read a test matrix
if (ReadMatrixMarket("test-symmetric.mtx", m)):
	print "Matrix read successfully!"
else:
	print "Error reading matrix!"

# Print the contents of m
print m

if (WriteMatrixMarket("test-symmetric-out.mtx", m, True)):
	print "Matrix written successfully!"
else:
	print "Error writing matrix!"

