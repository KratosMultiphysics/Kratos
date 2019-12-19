from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics import IsDistributedRun

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.RANSApplication.trilinos_fluid_solver_no_replace import TrilinosFluidSolverNoReplace as fluid_solver_no_replace
elif (not IsDistributedRun()):
    from KratosMultiphysics.RANSApplication.fluid_solver_no_replace import FluidSolverNoReplace as fluid_solver_no_replace
else:
    raise Exception("Distributed run requires TrilinosApplication")


class PeriodicFluidDynamicsAnalysis(FluidDynamicsAnalysis):
    def _CreateSolver(self):
        return fluid_solver_no_replace(self.model, self.project_parameters["solver_settings"])
