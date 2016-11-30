from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
import time
CheckForPreviousImport()




class RBFMapper:
    def __init__(self, master_model_part, slave_model_part, radius, order=2):
        self.masterModelPart = master_model_part
        self.slaveModelPart  = slave_model_part
        self.order = order
        self.radius = radius#
        # definition of the solvers
        tol = 1e-6
        max_it = 1000
        verbosity = 1
        m = 10
#        self.linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, tol, max_it, verbosity, m)
        self.linear_solver = SuperLUSolver()
 
        self.rbfProcess =  CustomRbfMapperProcess(self.linear_solver, self.masterModelPart, self.slaveModelPart, self.radius, self.order)
        
    def MapFromMasterToSlave(self,masterField, slaveField):
        print("Mapping from MASTER to SLAVE ...............")
        status = self.rbfProcess.MapFromMasterToSlave(masterField, slaveField)
        print("Mapping from MASTER to SLAVE Completed", status)
        time.sleep(2)
        return status
        
    def MapFromSlaveToMaster(self,slaveField, masterField):
        status = self.rbfProcess.MapFromSlaveToMaster(slaveField, masterField)
        print("Mapping from SLAVE to MASTER Completed")
        return status
        
