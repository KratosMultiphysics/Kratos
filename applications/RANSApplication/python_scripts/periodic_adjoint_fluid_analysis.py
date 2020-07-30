from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import IsDistributedRun

from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis
if (IsDistributedRun()):
    raise Exception("Distributed runs are not yet supported with periodic adjoint analysis")
else:
    from KratosMultiphysics.RANSApplication.adjoint_fluid_solver_no_replace import AdjointFluidSolverNoReplace as adjoint_fluid_solver_no_replace


class PeriodicAdjointFluidAnalysis(AdjointFluidAnalysis):
    def _CreateSolver(self):
        return adjoint_fluid_solver_no_replace(self.model, self.project_parameters["solver_settings"])
