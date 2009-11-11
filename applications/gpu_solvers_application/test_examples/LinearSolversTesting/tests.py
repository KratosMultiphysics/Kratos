import time
import fluid_only_var

#including kratos path
kratos_libs_path            = fluid_only_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = fluid_only_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosGPUSolversApplication import *

space_utils = UblasSparseSpace()

pA  = space_utils.CreateEmptyMatrixPointer()
pb = space_utils.CreateEmptyVectorPointer()

A =  (pA).GetReference()
b = (pb).GetReference()




#matrix_filename = "mat8.mm"
matrix_filename = "A_0.01005_.mm"
ReadMatrixMarketMatrix(matrix_filename,A)
print "finished reading matrix"

#vector_filename = "vecb8.mm"
vector_filename = "b_0.01005_.mm"
ReadMatrixMarketVector(vector_filename,b)
print "finished reading vector"

print "reading finished"

preSweeps = Vector(5)
postSweeps = Vector(5)

preSweeps[0] = 1
postSweeps[0] = 1
preSweeps[1] = postSweeps[1] = 1
preSweeps[2] = postSweeps[2] = 1
preSweeps[3] = postSweeps[3] = 1
preSweeps[4] = postSweeps[4] = 1

W = 4.0/3.0

#TESTING IF CG AND BICGSTAB CAN ITERATE INDEFINITIVELY
x = Vector(len(b))
for i in range(1, 2):
	print "----- iteration ", i, " ------"
	space_utils.SetToZeroVector(x)	
	
	t1 = time.time()
	#precond = KratosAMGPreconditioner(W, 2, True, 5, 1000, preSweeps, postSweeps, True)
	#precond = DiagonalpreconditionerType()
	linear_solver = GPUCGSolver(1e-9, 5000)#, precond)
	linear_solver.Solve(A, x, b)
	print linear_solver
	print "\n\n"
	t2 = time.time()
	print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
	linear_solver=0
	

exit(0)

