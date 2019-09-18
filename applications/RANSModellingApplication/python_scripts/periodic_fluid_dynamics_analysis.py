from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.RANSModellingApplication.fluid_solver_no_replace import FluidSolverNoReplace
from KratosMultiphysics.RANSModellingApplication.trilinos_fluid_solver_no_replace import TrilinosFluidSolverNoReplace

class PeriodicFluidDynamicsAnalysis(FluidDynamicsAnalysis):
    def _CreateSolver(self):
        if self.project_parameters["problem_data"]["parallel_type"].GetString() == "OpenMP":
            return FluidSolverNoReplace(self.model,self.project_parameters["solver_settings"])
        else:
            return TrilinosFluidSolverNoReplace(self.model,self.project_parameters["solver_settings"])

