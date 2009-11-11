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


fileOutput = open('output_preconditioner_test.log', 'w')


#matrix_filename = "mat8.mm"
matrix_filename = "A_0.01005_.mm"
ReadMatrixMarketMatrix(matrix_filename,A)
print "finished reading matrix"

#vector_filename = "vecb8.mm"
vector_filename = "b_0.01005_.mm"
ReadMatrixMarketVector(vector_filename,b)
print "finished reading vector"

print "reading finished"


maxSweeps = 4
maxLevelOfHierarchy = 5
minSizeOfLastLevel = 1000
maxLevelsRoh = 4
W = 4.0 / 3.0

preSweeps = Vector(maxLevelOfHierarchy)
postSweeps = Vector(maxLevelOfHierarchy)
preSweeps[0] = 1
postSweeps[0] = 1
preSweeps[1] = 2
postSweeps[1] = 2
preSweeps[2] = 2
postSweeps[2] = 2
preSweeps[3] = 2
postSweeps[3] = 2
preSweeps[4] = 2
postSweeps[4] = 2

roh = 2

#print 'GPU CG SOLVER\n'
#x = Vector(len(b))
#space_utils.SetToZeroVector(x)
#inst1 = time.time()
#linear_solver = GPUCGSolver(1e-9, 5000)
#inst2 = time.time()
#linear_solver.Solve(A, x, b)
#inst3 = time.time()
#print linear_solver
#print 'Creation time: %0.3f ms' % ((inst2-inst1)*1000.0)
#print 'Solve time: %0.3f ms' % ((inst3-inst2)*1000.0)
#linear_solver = 0
#print '\n\n'

for iter in range(1, 2):
	print ' --- ITERATION %i --- \n' %iter
	inst1 = time.time()
	#precond = KratosAMGPreconditioner(W, roh, True, maxLevelOfHierarchy, minSizeOfLastLevel, preSweeps, postSweeps)
	x = Vector(len(b))
	space_utils.SetToZeroVector(x)
	#linear_solver = GPUCGSolver(1e-9, 5000, precond)
	linear_solver = AMGSolver(1e-9, 5000, W, roh, False, maxLevelOfHierarchy, minSizeOfLastLevel, preSweeps, postSweeps)
	inst2 = time.time()											
	linear_solver.Solve(A, x, b)
	inst3 = time.time()
	print linear_solver
	print 'Creation time: %0.3f ms' % ((inst2-inst1)*1000.0)
	print 'Solve time: %0.3f ms' % ((inst3-inst2)*1000.0)
	linear_solver = 0

exit(0)

fileOutput.write('	Showing solving parameters and time for GPUCG + preconditioner\n\n')
fileOutput.write('LvlRoh | preSweep1 | preSweep2 | preSweep3 | postSweep1 | postSweep2 | postSweep3 | creationTime | solvingTime\n')

for tol in [1e-9, 1e-6, 1e-3]:
	for roh in range(maxLevelsRoh, 0, -1):

		preSweeps = Vector(maxLevelOfHierarchy)
		postSweeps = Vector(maxLevelOfHierarchy)
	
		#first level		
		for pre_sweeps1 in range(maxSweeps, 0, -1):
			preSweeps[0] = postSweeps[0] = pre_sweeps1
			#second level
			for pre_sweeps2 in range(maxSweeps, 0, -1):
				preSweeps[1] = postSweeps[1] = pre_sweeps2
				#third level
#				for pre_sweeps3 in range(maxSweeps, 0, -1):
#					preSweeps[2] = postSweeps[2] = pre_sweeps3						
				print "\n		---New Test---"
				print "LvlRoh ", roh, ", preSweep1 ", preSweeps[0], ", postSweep1 ", postSweeps[0], ", preSweep2 ", preSweeps[1], ", postSweep2 ", postSweeps[1]#, ", preSweep3 ", preSweeps[2], ", postSweep3 ", postSweeps[2], ", preSweep4 ", pre_sweeps4, ", postSweep4 ", post_sweeps4, ", preSweep5 ", pre_sweeps5, ", postSweep5 ", post_sweeps5
				inst1 = time.time()
				precond = KratosAMGPreconditioner(W, roh, True, maxLevelOfHierarchy, minSizeOfLastLevel, preSweeps, postSweeps)
				x = Vector(len(b))
				space_utils.SetToZeroVector(x)
				linear_solver = GPUCGSolver(tol, 5000, precond)
				inst2 = time.time()											
				linear_solver.Solve(A, x, b)
				inst3 = time.time()
				print linear_solver
				print 'Creation time: %0.3f ms' % ((inst2-inst1)*1000.0)
				print 'Solve time: %0.3f ms' % ((inst3-inst2)*1000.0)
				linear_solver = 0
				precond = 0
				fileOutput.write('%i | %i | %i | %i | %i | %f | %f\n' %(roh ,preSweeps[0] ,preSweeps[1] ,postSweeps[0] ,postSweeps[1] ,((inst2-inst1)*1000.0) ,((inst3-inst2)*1000.0)) )
							




fileOutput.close()
#exit(0)

print "		---GPUCG test---"
inst1 = time.time()
linear_solver =  GPUCGSolver(1e-9, 5000)
x = Vector(len(b))
space_utils.SetToZeroVector(x)
inst2 = time.time()
linear_solver.Solve(A,x,b)
inst3 = time.time()
print linear_solver
print 'Creation time: %0.3f ms' % ((inst2-inst1)*1000.0)
print 'Solve time: %0.3f ms' % ((inst3-inst2)*1000.0)

print "		---CG test---"
inst1 = time.time()
linear_solver =  CGSolver(1e-9, 5000)
x = Vector(len(b))
space_utils.SetToZeroVector(x)
inst2 = time.time()
linear_solver.Solve(A,x,b)
inst3 = time.time()
print linear_solver
print 'Creation time: %0.3f ms' % ((inst2-inst1)*1000.0)
print 'Solve time: %0.3f ms' % ((inst3-inst2)*1000.0)



