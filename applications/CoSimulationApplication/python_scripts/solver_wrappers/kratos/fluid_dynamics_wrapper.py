from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from . import kratos_base_wrapper

# Other imports
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

def CreateSolver(model, settings, solver_name):
    return KratosFluidSolver(model, settings, solver_name)

class KratosFluidSolver(kratos_base_wrapper.KratosBaseInterface):
    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysis(self.model, self.project_parameters)
