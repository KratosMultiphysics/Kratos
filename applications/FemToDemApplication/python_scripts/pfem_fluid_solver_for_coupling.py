from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_solver as pfem_fluid_solver
import KratosMultiphysics.FemToDemApplication.pfem_check_and_prepare_model_process_fluid_for_coupling as pfem_check_and_prepare_model_process_fluid_for_coupling

def CreateSolver(model, FEM_model_part, parameters):
    return PfemFluidSolverForCoupling(model, FEM_model_part, parameters)

#============================================================================================================================
class PfemFluidSolverForCoupling(pfem_fluid_solver.PfemFluidSolver):
#============================================================================================================================
    def __init__(self, model, FEM_model_part, parameters):
        self.FEM_model_part = FEM_model_part
        super(PfemFluidSolverForCoupling,self).__init__(model, parameters)

#============================================================================================================================
    def CheckAndPrepareModelProcess(self, params):
        # CheckAndPrepareModelProcess creates the fluid_computational model part
        pfem_check_and_prepare_model_process_fluid_for_coupling.CheckAndPrepareModelProcessForCoupling(self.main_model_part,
                                                                                                       self.FEM_model_part,
                                                                                                       params).Execute()
#============================================================================================================================