import time

from KratosMultiphysics import *
from KratosMultiphysics.GPUSolversApplication import *

space_utils = UblasSparseSpace()

pA  = space_utils.CreateEmptyMatrixPointer()
pb = space_utils.CreateEmptyVectorPointer()

A =  (pA).GetReference()
b = (pb).GetReference()




matrix_filename = "mat8.mm"
#matrix_filename = "A_0.01005_.mm"
#matrix_filename = "toTest/mat80.mm"

ReadMatrixMarketMatrix(matrix_filename,A)
print "finished reading matrix"

vector_filename = "vecb8.mm"
#vector_filename = "b_0.01005_.mm"
#vector_filename = "toTest/vecb80.mm"

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

space_utils.SetToZeroVector(x)	

print "AMG standalone"
t1 = time.time()
linear_solver2 = AMGSolver(1e-9, 5000, W, 2, False, 5, 100, preSweeps, postSweeps)
linear_solver2.Solve(A, x, b)
print linear_solver2
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)

exit(0)

space_utils.SetToZeroVector(x)	

print "GPU CG + SA-AMG"
t1 = time.time()
precond1 = KratosAMGPreconditioner(W, 2, True, 5, 1000, preSweeps, postSweeps, True)
linear_solver1 = GPUCGSolver(1e-9, 5000, precond1)
linear_solver1.Solve(A, x, b)
print linear_solver1
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver=0


space_utils.SetToZeroVector(x)	

print "GPU CG + diag precond"
t1 = time.time()
precond2 = GPUDiagonalPreconditioner()
linear_solver2 = GPUCGSolver(1e-9, 5000, precond2)
linear_solver2.Solve(A, x, b)
print linear_solver2
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver=0

space_utils.SetToZeroVector(x)	

print "GPU CG"
t1 = time.time()
linear_solver2 = GPUCGSolver(1e-9, 5000)
linear_solver2.Solve(A, x, b)
print linear_solver2
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver=0


space_utils.SetToZeroVector(x)	

print "CPU CG + diag precond"
t1 = time.time()
precond3 = DiagonalPreconditioner()
linear_solver3 = CGSolver(1e-9, 5000, precond3)
linear_solver3.Solve(A, x, b)
print linear_solver3
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver3=0

space_utils.SetToZeroVector(x)	

print "GPU BICG + diag precond"
t1 = time.time()
precond2 = GPUDiagonalPreconditioner()
linear_solver2 = GPUBICGSTABSolver(1e-9, 5000, precond2)
linear_solver2.Solve(A, x, b)
print linear_solver2
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver=0

space_utils.SetToZeroVector(x)	

print "GPU BICG"
t1 = time.time()
linear_solver2 = GPUBICGSTABSolver(1e-9, 5000)
linear_solver2.Solve(A, x, b)
print linear_solver2
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver=0

space_utils.SetToZeroVector(x)	

print "CPU BICG + diag precond"
t1 = time.time()
precond2 = DiagonalPreconditioner()
linear_solver2 = BICGSTABSolver(1e-9, 5000, precond2)
linear_solver2.Solve(A, x, b)
print linear_solver2
print "\n\n"
t2 = time.time()
print 'Solve time: %0.3f ms' % ((t2-t1)*1000.0)
#linear_solver=0


#linear_solver=0

exit(0)

